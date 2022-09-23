/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2010,2011,2012,2013 TELEMATICS LAB, Politecnico di Bari
 *
 * This file is part of LTE-Sim
 *
 * LTE-Sim is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as
 * published by the Free Software Foundation;
 *
 * LTE-Sim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LTE-Sim; if not, see <http://www.gnu.org/licenses/>.
 *
 * Author: Giuseppe Piro <g.piro@poliba.it>
 *         Lukasz Rajewski <lukasz.rajewski@gmail.com> (optimized PRB allocation)
 */


#include "downlink-nvs-scheduler.h"
#include "../mac-entity.h"
#include "../../packet/Packet.h"
#include "../../packet/packet-burst.h"
#include "../../../device/NetworkNode.h"
#include "../../../flows/radio-bearer.h"
#include "../../../protocolStack/rrc/rrc-entity.h"
#include "../../../flows/application/Application.h"
#include "../../../device/ENodeB.h"
#include "../../../protocolStack/mac/AMCModule.h"
#include "../../../phy/lte-phy.h"
#include "../../../core/spectrum/bandwidth-manager.h"
#include "../../../flows/MacQueue.h"
#include "../../../utility/eesm-effective-sinr.h"
#include <cstdio>
#include <limits>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <cstdlib>

DownlinkNVSScheduler::DownlinkNVSScheduler(std::string config_fname, bool is_optimal)
  : is_optimal_(is_optimal) {
  std::ifstream ifs(config_fname, std::ifstream::in);
  if (ifs.is_open()) {
    std::string line;
    while (std::getline(ifs, line)) {
      if (line[0] == '#') {
        continue;
      }
      std::string args;
      std::istringstream iss(line);
      iss >> args;
      if (args == "num_slices:") {
        iss >> num_slices_;
      }
      else if (args == "slice_weights:") {
        for (int i = 0; i < num_slices_; ++i) {
          iss >> slice_weights_[i];
        }
      }
      else if (args == "slice_ues:") {
        int begin_id = 0, num_ue;
        for (int i = 0; i < num_slices_; ++i) {
          iss >> num_ue;
          for (int j = 0; j < num_ue; ++j) {
            user_to_slice_[begin_id + j] = i;
          }
          begin_id += num_ue;
        }
      }
      else if (args == "slice_algos:") {
        int algo;
        for (int i = 0; i < num_slices_; ++i) {
          iss >> algo;
          if (algo == 0)
            slice_algo_[i] = MT;
          else if (algo == 1)
            slice_algo_[i] = PF;
          else if (algo == 2)
            slice_algo_[i] = MLWDF;
          else
            slice_algo_[i] = PF;
        }
      }
    }
  }
  else {
    throw std::runtime_error("Fail to open configuration file.");
  }
  ifs.close();

  SetMacEntity (0);
  CreateUsersToSchedule();
}

DownlinkNVSScheduler::~DownlinkNVSScheduler()
{
  Destroy ();
}

int DownlinkNVSScheduler::SelectSliceToServe()
{
  int slice_id = 0;
  double max_score = 0;

#ifdef SCHEDULER_DEBUG
  for(int i = 0; i < num_type2_slices_; ++i) {
    std::cout << i << ": weight: " << type2_weights_[i] << " time: " << type2_exp_time_[i] << std::endl;
  }
#endif

  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();
  std::vector<bool> slice_with_queue(num_slices_, false);
  for (auto it = bearers->begin(); it != bearers->end(); it++) {
    RadioBearer *bearer = (*it);
    int user_id = bearer->GetUserID();

    if (bearer->HasPackets () && bearer->GetDestination ()->GetNodeState () == NetworkNode::STATE_ACTIVE)
		{
		  int dataToTransmit;
		  if (bearer->GetApplication ()->GetApplicationType () == Application::APPLICATION_TYPE_INFINITE_BUFFER) {
			  dataToTransmit = 100000000;
			}
		  else {
			  dataToTransmit = bearer->GetQueueSize ();
			}
      if (dataToTransmit > 0) {
        slice_with_queue[user_to_slice_[user_id]] = true;
      }
    }
  }
  for (int i = 0; i < num_slices_; ++i) {
    if (! slice_with_queue[i]) continue;
    if (slice_exp_time_[i] == 0) {
      slice_id = i;
      break;
    }
    else {
      double score = slice_weights_[i] / slice_exp_time_[i];
      if (score >= max_score) {
        max_score = score;
        slice_id = i;
      }
    }
  }
  for (int i = 0; i < num_slices_; ++i) {
    if (! slice_with_queue[i]) continue;
    slice_exp_time_[i] = (1-beta_) * slice_exp_time_[i];
    if (i == slice_id) {
      slice_exp_time_[i] += beta_ * 1;
    }
  }
  return slice_id;
}

void DownlinkNVSScheduler::SelectFlowsToSchedule (int slice_serve)
{
  ClearUsersToSchedule();

  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  memset(slice_priority_, 0, sizeof(int) * MAX_SLICES);
  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
	{
	  //SELECT FLOWS TO SCHEDULE
	  RadioBearer *bearer = (*it);
    int user_id = bearer->GetUserID();

    if (user_to_slice_[user_id] != slice_serve)
      continue;

	  if (bearer->HasPackets () && bearer->GetDestination ()->GetNodeState () == NetworkNode::STATE_ACTIVE) {
		  //compute data to transmit
		  int dataToTransmit;
		  if (bearer->GetApplication ()->GetApplicationType () == Application::APPLICATION_TYPE_INFINITE_BUFFER) {
			  dataToTransmit = 100000000;
			}
		  else {
			  dataToTransmit = bearer->GetQueueSize ();
			}

		  //compute spectral efficiency
		  ENodeB *enb = (ENodeB*) GetMacEntity ()->GetDevice ();
		  ENodeB::UserEquipmentRecord *ueRecord = enb->GetUserEquipmentRecord (bearer->GetDestination ()->GetIDNetworkNode ());
		  std::vector<double> spectralEfficiency;
		  std::vector<int> cqiFeedbacks = ueRecord->GetCQI ();
		  int numberOfCqi = cqiFeedbacks.size ();
		  for (int i = 0; i < numberOfCqi; i++) {
			  double sEff = GetMacEntity ()->GetAmcModule ()->GetEfficiencyFromCQI (cqiFeedbacks.at (i));
			  spectralEfficiency.push_back (sEff);
			}

		  //create flow to scheduler record
      int slice_id = user_to_slice_[user_id];
      if (bearer->GetPriority() > slice_priority_[slice_id]) {
        slice_priority_[slice_id] = bearer->GetPriority();
      }
      InsertFlowToUser(bearer, dataToTransmit, spectralEfficiency, cqiFeedbacks);
    }
	}
}

void
DownlinkNVSScheduler::DoSchedule (void)
{
#ifdef SCHEDULER_DEBUG
	std::cout << "\nStart DL packet scheduler for node "
			<< GetMacEntity ()->GetDevice ()->GetIDNetworkNode()
      << " ts: " << GetTimeStamp() << std::endl;
#endif

  // some logic to determine the slice to serve
  int slice_serve = SelectSliceToServe();
  UpdateAverageTransmissionRate (slice_serve);
  SelectFlowsToSchedule (slice_serve);

  if (GetUsersToSchedule()->size() != 0) {
    if (is_optimal_)
      RBsAllocationOptimalPF();
    else
      RBsAllocation();
  }

  StopSchedule ();
}

void
DownlinkNVSScheduler::DoStopSchedule(void)
{
  PacketBurst* pb = new PacketBurst();
  UsersToSchedule *uesToSchedule = GetUsersToSchedule();
  for (auto it = uesToSchedule->begin(); it != uesToSchedule->end(); it++) {
    UserToSchedule* user = *it;
    int availableBytes = user->GetAllocatedBits() / 8;
    // let's not reallocate RBs between users firstly
    // when the flow of highest priority has no data, the left availabe bytes
    // are reallocated to lower prioritized flows
    for (int i = MAX_BEARERS-1; i >= 0; i--) {
      if (availableBytes <= 0)
        break;
      if (user->m_dataToTransmit[i] > 0) {
        assert(user->m_bearers[i] != NULL);
        int dataTransmitted = min(availableBytes, user->m_dataToTransmit[i]);
        availableBytes -= dataTransmitted;
        user->m_bearers[i]->UpdateTransmittedBytes(
            dataTransmitted
            );
        user->m_bearers[i]->UpdateCumulateRBs(
            user->GetListOfAllocatedRBs()->size()
            );
        std::cerr << GetTimeStamp()
          << " app: " << user->m_bearers[i]->GetApplication()->GetApplicationID()
          << " bytes: " << user->m_bearers[i]->GetCumulateBytes()
          << " rbs: " << user->m_bearers[i]->GetCumulateRBs()
          << " delay: " << user->m_bearers[i]->GetHeadOfLinePacketDelay()
          << " user: " << user->GetUserID()
          << " slice: " << user_to_slice_[user->GetUserID()]
          << std::endl;

        RlcEntity *rlc = user->m_bearers[i]->GetRlcEntity ();
        PacketBurst* pb2 = rlc->TransmissionProcedure (dataTransmitted);

        if (pb2->GetNPackets () > 0)
        {
          std::list<Packet*> packets = pb2->GetPackets ();
          std::list<Packet* >::iterator it;
          for (it = packets.begin (); it != packets.end (); it++)
          {
            Packet *p = (*it);
            pb->AddPacket (p->Copy ());
          }
        }
        delete pb2;  
      }
    }
  }
  UpdateTimeStamp();

  GetMacEntity ()->GetDevice ()->SendPacketBurst (pb);
}

void
DownlinkNVSScheduler::RBsAllocation()
{
  int nb_rbs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  size_t rbg_size = get_rbg_size(nb_rbs);
  // currently nb_rbgs should be divisible
  nb_rbs = nb_rbs - (nb_rbs % rbg_size);
  int nb_rbgs = nb_rbs / rbg_size;
  UsersToSchedule* users = GetUsersToSchedule();
  double metrics[nb_rbgs][users->size()];
  for (size_t i = 0; i < nb_rbgs; i++) {
    for (size_t j = 0; j < users->size(); j++) {
      metrics[i][j] = ComputeSchedulingMetric(
        users->at(j),
        users->at(j)->GetSpectralEfficiency().at(i * rbg_size)
      );
    }
  }
  for (int i = 0; i < nb_rbgs; i++) {
    double targetMetric = std::numeric_limits<double>::lowest();
    UserToSchedule * scheduledUE = NULL;
    for (size_t j = 0; j < users->size(); ++j) {
      UserToSchedule* user = users->at(j);
      // allocate the rbg if the user has highest metric and data to transmit
      if (metrics[i][j] > targetMetric &&
          user->GetListOfAllocatedRBs()->size() < user->m_requiredRBs ) {
        targetMetric = metrics[i][j];
        scheduledUE = user;
      }
    }
    if (scheduledUE) {
      int l = i*rbg_size, r = (i+1)*rbg_size;
      for (int i = l; i < r; ++i) {
        scheduledUE->GetListOfAllocatedRBs()->push_back(i);
      }
    }
  }
  AMCModule *amc = GetMacEntity()->GetAmcModule();
  PdcchMapIdealControlMessage *pdcchMsg = new PdcchMapIdealControlMessage();
  std::cout << GetTimeStamp() << std::endl;
  for (size_t j = 0; j < users->size(); j++) {
    UserToSchedule* ue = users->at(j);

    if (ue->GetListOfAllocatedRBs()->size() > 0) {
      std::vector<double> estimatedSinrValues;

      std::cout << "User(" << ue->GetUserID() << ") allocated RBGS:";
      for (size_t i = 0; i < ue->GetListOfAllocatedRBs()->size (); i++ ) {
        int rbid = ue->GetListOfAllocatedRBs()->at(i);
        if (rbid % rbg_size == 0)
          std::cout << " " << rbid / rbg_size <<
            "(" << ue->GetCqiFeedbacks().at(rbid) << ")";
        
        double sinr = amc->GetSinrFromCQI(
          ue->GetCqiFeedbacks().at(
          ue->GetListOfAllocatedRBs()->at(i)));
        estimatedSinrValues.push_back (sinr);
      }
      double effectiveSinr = GetEesmEffectiveSinr(estimatedSinrValues);
      std::cout << " final_cqi: " << amc->GetCQIFromSinr(effectiveSinr) << std::endl;
      int mcs = amc->GetMCSFromCQI(amc->GetCQIFromSinr(effectiveSinr));
      int transportBlockSize = amc->GetTBSizeFromMCS(mcs, ue->GetListOfAllocatedRBs()->size());

      // // calculate the transport block size as if the user can have multiple mcs
      //int mcs = 1;
      //int transportBlockSize = 0;
      //for (int i = 0; i < estimatedSinrValues.size(); i++) {
      //    transportBlockSize += amc->GetTBSizeFromMCS(amc->GetMCSFromCQI(amc->GetCQIFromSinr(estimatedSinrValues[i])), 1);
      //}
      ue->UpdateAllocatedBits(transportBlockSize);
      for (size_t rb = 0; rb < ue->GetListOfAllocatedRBs()->size(); rb++) {
        pdcchMsg->AddNewRecord(
            PdcchMapIdealControlMessage::DOWNLINK,
            ue->GetListOfAllocatedRBs()->at(rb),
            ue->GetUserNode(),
            mcs
          );
      }
    }
  }
  if (pdcchMsg->GetMessage()->size() > 0) {
    GetMacEntity()->GetDevice()->GetPhy()->SendIdealControlMessage(pdcchMsg);
  }
  delete pdcchMsg;
}

double
DownlinkNVSScheduler::ComputeSchedulingMetric(UserToSchedule* user, double spectralEfficiency)
{
  double metric = 0;
  int slice_id = user_to_slice_[user->GetUserID()];
  // get the aggregate transmission rate of all flows
  double average_rate = user->GetAverageTransmissionRate();

  if (schedule_scheme_ == 0) {
    switch (slice_algo_[slice_id]) {
      case MT:
        metric = spectralEfficiency;
        break;
      case PF:
        metric = (spectralEfficiency * 180000.) / average_rate;
        break;
      default:
        metric = spectralEfficiency;
    }
  }
  else if (schedule_scheme_ == 1) {
    // select the highest priority first
    // if the user doesn't have data with slice_priority_, it's not scheduled at all
    if (user->m_dataToTransmit[slice_priority_[slice_id]] == 0){
      metric = 0;
    }
    else {
      if (slice_algo_[slice_id] == MT) {
        metric = spectralEfficiency;
      }
      else if (slice_algo_[slice_id] == PF) {
        metric = (spectralEfficiency * 180000.) / average_rate;
      }
      else if (slice_algo_[slice_id] == MLWDF) {
        RadioBearer* bearer = user->m_bearers[slice_priority_[slice_id]];
        double HoL = bearer->GetHeadOfLinePacketDelay();
        metric = HoL * (spectralEfficiency * 180000.) / average_rate;
      }
      else {
        throw std::runtime_error("Invalid Scheduling algorithm");
      }
    }
  }
  return metric;
}

void
DownlinkNVSScheduler::UpdateAverageTransmissionRate (int slice_serve)
{
  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
  {
    // RadioBearer *bearer = (*it);
    // int user_id = bearer->GetUserID();
    // if (user_to_slice_[user_id] != slice_serve)
    //   continue;

	  RadioBearer *bearer = (*it);
    bearer->UpdateAverageTransmissionRate ();
  }
}

void
DownlinkNVSScheduler::RBsAllocationOptimalPF() {
  int nb_rbs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  size_t rbg_size = get_rbg_size(nb_rbs);
  // currently nb_rbgs should be divisible
  nb_rbs = nb_rbs - (nb_rbs % rbg_size);
  int nb_rbgs = nb_rbs / rbg_size;

  UsersToSchedule* users = GetUsersToSchedule();
  std::vector<int> user_highest_cqi;
  for (auto it = users->begin(); it != users->end(); it++) {
    int highest_cqi = 0;
    for (int i = 0; i < nb_rbgs; i++) {
      int cqi = (*it)->GetCqiFeedbacks().at(i * rbg_size);
      if (highest_cqi < cqi)
        highest_cqi = cqi;
    }
    user_highest_cqi.push_back(highest_cqi);
  }
  std::vector<int> best_assigned_mcs;
  std::vector<UserToSchedule*> best_rbgs_assignment;
  double highest_pf_metric = 0;
  int cqi_search_range = 4;
  int num_sample = 300;
  for (int i = 0; i < num_sample; i++) {
    std::vector<int> assigned_mcs;
    std::vector<UserToSchedule*> rbgs_assignment(nb_rbgs, NULL);
    double pf_metric = 0;
    // std::cout << "MCS(highest_cqi): ";
    for (size_t i = 0; i < user_highest_cqi.size(); i++) {
      assigned_mcs.push_back(
        max( user_highest_cqi[i] - rand() % cqi_search_range, 1 )
      );
      // std::cout << assigned_mcs.back() << "(" << user_highest_cqi[i] << ") ";
    }
    // std::cout << std::endl;
    pf_metric = AssignRBsGivenMCS(assigned_mcs, rbgs_assignment);
    if (highest_pf_metric < pf_metric) {
      highest_pf_metric = pf_metric;
      best_assigned_mcs = assigned_mcs;
      best_rbgs_assignment = rbgs_assignment;
    }
  }
  for (size_t rbg_id = 0; rbg_id < best_rbgs_assignment.size(); rbg_id++) {
    UserToSchedule* user = best_rbgs_assignment[rbg_id];
    assert(user != NULL);
    int l = rbg_id * rbg_size, r = (rbg_id+1) * rbg_size;
    for (int j = l; j < r; ++j) {
      user->GetListOfAllocatedRBs()->push_back(j);
    }
  }
  AMCModule *amc = GetMacEntity()->GetAmcModule();
  PdcchMapIdealControlMessage *pdcchMsg = new PdcchMapIdealControlMessage();

  std::cout << GetTimeStamp() << std::endl;
  for (auto it = users->begin(); it != users->end(); it++) {
    UserToSchedule* ue = *it;
    if (ue->GetListOfAllocatedRBs()->size() > 0) {
      std::vector<double> estimatedSinrValues;
      std::cout << "User(" << ue->GetUserID() << ") allocated RBGS:";

      for (size_t i = 0; i < ue->GetListOfAllocatedRBs()->size (); i++ ) {
        int rbid = ue->GetListOfAllocatedRBs()->at(i);
        if (rbid % rbg_size == 0)
          std::cout << " " << rbid / rbg_size <<
            "(" << ue->GetCqiFeedbacks().at(rbid) << ")";
        double sinr = amc->GetSinrFromCQI ( ue->GetCqiFeedbacks().at (rbid) );
        estimatedSinrValues.push_back (sinr);
      }
      double effectiveSinr = GetEesmEffectiveSinr(estimatedSinrValues);

      std::cout << " final_cqi: " << amc->GetCQIFromSinr(effectiveSinr) << std::endl;

      int mcs = amc->GetMCSFromCQI(amc->GetCQIFromSinr(effectiveSinr));
      int transportBlockSize = amc->GetTBSizeFromMCS(mcs, ue->GetListOfAllocatedRBs()->size());

      ue->UpdateAllocatedBits(transportBlockSize);
      for (size_t rb = 0; rb < ue->GetListOfAllocatedRBs()->size(); rb++) {
        pdcchMsg->AddNewRecord(
            PdcchMapIdealControlMessage::DOWNLINK,
            ue->GetListOfAllocatedRBs()->at(rb),
            ue->GetUserNode(),
            mcs
          );
      }
    }
  }
  if (pdcchMsg->GetMessage()->size() > 0) {
    GetMacEntity()->GetDevice()->GetPhy()->SendIdealControlMessage(pdcchMsg);
  }
  delete pdcchMsg;
}

double
DownlinkNVSScheduler::AssignRBsGivenMCS(std::vector<int>& assigned_mcs,
    std::vector<UserToSchedule*>& rbgs_assignment) {
  int nb_rbs = GetMacEntity()->GetDevice()->GetPhy() \
      ->GetBandwidthManager()->GetDlSubChannels().size();
  size_t rbg_size = get_rbg_size(nb_rbs);
  int nb_rbgs = rbgs_assignment.size();

  UsersToSchedule* users = GetUsersToSchedule();
  double pf_metric = 0;

  for (int i = 0; i < nb_rbgs; i++) {
    double highest_metric = -1;
    for (auto it = users->begin(); it != users->end(); it++) {
      int index = std::distance(users->begin(), it);
      int cqi = (*it)->GetCqiFeedbacks().at(i * rbg_size);
      double metric = 0;
      // the metric is not zero if cqi is gte the mcs
      if (assigned_mcs[index] <= cqi) {
        double sEff = GetMacEntity()->GetAmcModule()->GetEfficiencyFromCQI (assigned_mcs[index]);
        metric = sEff * 180000 / (*it)->GetAverageTransmissionRate();
      }
      if (highest_metric < metric) {
        highest_metric = metric;
        rbgs_assignment[i] = *it;
      }
    }
    pf_metric += highest_metric;
  }
  return pf_metric;
}
