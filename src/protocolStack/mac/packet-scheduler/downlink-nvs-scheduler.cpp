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
#include <cassert>
#include <cstring>

DownlinkNVSScheduler::DownlinkNVSScheduler(std::string config_fname)
{
  std::ifstream ifs(config_fname, std::ifstream::in);
  if (ifs.is_open()) {
    int begin_id = 0;
    ifs >> num_slices_;

    for (int i = 0; i < num_slices_; ++i)
      ifs >> slice_weights_[i];
    for (int i = 0; i < num_slices_; ++i) {
      slice_exp_time_[i] = slice_weights_[i];
      int num_ue;
      ifs >> num_ue;
      for (int j = 0; j < num_ue; ++j) {
        user_to_slice_[begin_id + j] = i;
      }
      begin_id += num_ue;
    }
    for (int i = 0; i < num_slices_; ++i) {
      int algo;
      ifs >> algo;
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
  else {
    throw std::runtime_error("Fail to open configuration file");
  }
  ifs.close();

  SetMacEntity (0);
  CreateUsersToSchedule();
}

DownlinkNVSScheduler::~DownlinkNVSScheduler()
{
  Destroy ();
}

void DownlinkNVSScheduler::SelectSliceToServe(int& slice_id)
{
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

	  if (bearer->HasPackets () && bearer->GetDestination ()->GetNodeState () == NetworkNode::STATE_ACTIVE)
		{
		  //compute data to transmit
		  int dataToTransmit;
		  if (bearer->GetApplication ()->GetApplicationType () == Application::APPLICATION_TYPE_INFINITE_BUFFER)
			{
			  dataToTransmit = 100000000;
			}
		  else
			{
			  dataToTransmit = bearer->GetQueueSize ();
			}

		  //compute spectral efficiency
		  ENodeB *enb = (ENodeB*) GetMacEntity ()->GetDevice ();
		  ENodeB::UserEquipmentRecord *ueRecord = enb->GetUserEquipmentRecord (bearer->GetDestination ()->GetIDNetworkNode ());
		  std::vector<double> spectralEfficiency;
		  std::vector<int> cqiFeedbacks = ueRecord->GetCQI ();
		  int numberOfCqi = cqiFeedbacks.size ();
		  for (int i = 0; i < numberOfCqi; i++)
			{
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
  int slice_serve = 0;
  SelectSliceToServe(slice_serve);
  UpdateAverageTransmissionRate (slice_serve);
  SelectFlowsToSchedule (slice_serve);

  if (GetUsersToSchedule()->size() == 0) {
  }
  else {
	  RBsAllocation ();
    //RBsAllocationForUE();
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
    if (availableBytes > 0) {
      // let's not reallocate RBs between flows firstly
      for (int i = MAX_BEARERS-1; i >= 0; i--) {
        if (user->m_dataToTransmit[i] > 0) {
          assert(user->m_bearers[i] != NULL);
          user->m_bearers[i]->UpdateTransmittedBytes(
            min(availableBytes, user->m_dataToTransmit[i])
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
	        PacketBurst* pb2 = rlc->TransmissionProcedure (availableBytes);

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
  int nb_rbgs = (nb_rbs + rbg_size - 1) / rbg_size;
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
  for (size_t j = 0; j < users->size(); j++) {
    UserToSchedule* ue = users->at(j);

    if (ue->GetListOfAllocatedRBs()->size() > 0) {
      std::vector<double> estimatedSinrValues;
      for (size_t rb = 0; rb < ue->GetListOfAllocatedRBs()->size (); rb++ ) {
        double sinr = amc->GetSinrFromCQI (
        ue->GetCqiFeedbacks ().at (
          ue->GetListOfAllocatedRBs ()->at (rb)));
        estimatedSinrValues.push_back (sinr);
      }
      double effectiveSinr = GetEesmEffectiveSinr(estimatedSinrValues);
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

// double
// DownlinkNVSScheduler::ComputeSchedulingMetric(RadioBearer *bearer, double spectralEfficiency, int subChannel)
// {
//   double metric = 0;
//   switch (intra_sched_) {
//     case MT:
//       metric = spectralEfficiency;
//       break;
//     case PF:
//       metric = (spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate();
//       break;
//     case MLWDF:
//     {
//       double a = -log10 (0.01) / 1;
//       double HOL = bearer->GetHeadOfLinePacketDelay ();
//       metric = (a * HOL) * ((spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate ());
//       break;
//     }

//     default:
//       metric = spectralEfficiency;
//   }
//   return metric;
// }

double
DownlinkNVSScheduler::ComputeSchedulingMetric(UserToSchedule* user, double spectralEfficiency)
{
  double metric = 0;
  double averageRate = 1;
  int slice_id = user_to_slice_[user->GetUserID()];
  for (int i = 0; i < MAX_BEARERS; ++i) {
    if (user->m_bearers[i]) {
      averageRate += user->m_bearers[i]->GetAverageTransmissionRate();
    }
  }
  // recomputation of averageRate
  // double beta = 0.1;
  // AMCModule *amc = GetMacEntity()->GetAmcModule();
  // int bitsSent = amc->GetTBSizeFromMCS(
  //   amc->GetMCSFromCQI(user->GetWidebandCQI()),
  //   user->GetListOfAllocatedRBs()->size()); 
  // averageRate += (beta * bitsSent / 0.001);

  if (schedule_scheme_ == 0) {
    switch (slice_algo_[slice_id]) {
      case MT:
        metric = spectralEfficiency;
        break;
      case PF:
        metric = (spectralEfficiency * 180000.) / averageRate;
        break;
      default:
        metric = spectralEfficiency;
    }
  }
  else if (schedule_scheme_ == 1) {
    if (user->m_dataToTransmit[slice_priority_[slice_id]] == 0){
      metric = 0;
    }
    else {
      if (slice_algo_[slice_id] == MT) {
        metric = spectralEfficiency;
      }
      else if (slice_algo_[slice_id] == PF) {
        metric = (spectralEfficiency * 180000.) / averageRate;
      }
      else if (slice_algo_[slice_id] == MLWDF) {
        RadioBearer* bearer = user->m_bearers[slice_priority_[slice_id]];
        //assert(bearer->)
        double HoL = bearer->GetHeadOfLinePacketDelay();
        metric = HoL * (spectralEfficiency * 180000.) / averageRate;
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
	  RadioBearer *bearer = (*it);
    bearer->UpdateAverageTransmissionRate ();
  }
}
