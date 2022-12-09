/* project: RadioSaber; Mode: C++
 * Copyright (c) 2021, 2022, 2023, 2024 University of Illinois Urbana Champaign
 *
 * This file is part of RadioSaber, which is a project built upon LTE-Sim in 2022
 *
 * RadioSaber is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * RadioSaber is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author: Yongzhou Chen <yongzhouc@outlook.com>
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
#include "../../../load-parameters.h"
#include <jsoncpp/json/json.h>
#include <cstdio>
#include <limits>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <cstdlib>

DownlinkNVSScheduler::DownlinkNVSScheduler(std::string config_fname, bool is_nongreedy)
  : is_nongreedy_(is_nongreedy) {
  std::ifstream ifs(config_fname);
  if (!ifs.is_open()) {
    throw std::runtime_error("Fail to open configuration file.");
  }
  Json::Reader reader;
  Json::Value obj;
  reader.parse(ifs, obj);
  ifs.close();
  const Json::Value& ues_per_slice = obj["ues_per_slice"];
  num_slices_ = ues_per_slice.size();
  int num_ue;
  for (int i = 0; i < num_slices_; i++) {
    num_ue = ues_per_slice[i].asInt();
    for (int j = 0; j < num_ue; j++) {
      user_to_slice_.push_back(i);
    }
  }
  const Json::Value& slice_schemes = obj["slices"];
  for (int i = 0; i < slice_schemes.size(); i++) {
    int n_slices = slice_schemes[i]["n_slices"].asInt();
    for (int j = 0; j < n_slices; j++) {
      slice_weights_.push_back(
        slice_schemes[i]["weight"].asDouble()
      );
      slice_algo_params_.emplace_back(
        slice_schemes[i]["algo_alpha"].asInt(),
        slice_schemes[i]["algo_beta"].asInt(),
        slice_schemes[i]["algo_epsilon"].asInt(),
        slice_schemes[i]["algo_psi"].asInt()
      );
    }
  }
  slice_priority_.resize(num_slices_);
  std::fill(slice_priority_.begin(), slice_priority_.end(), 0);
  slice_ewma_time_.resize(num_slices_);
  std::fill(slice_ewma_time_.begin(), slice_ewma_time_.end(), 0);

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
    if (slice_ewma_time_[i] == 0) {
      slice_id = i;
      break;
    }
    else {
      double score = slice_weights_[i] / slice_ewma_time_[i];
      if (score >= max_score) {
        max_score = score;
        slice_id = i;
      }
    }
  }
  for (int i = 0; i < num_slices_; ++i) {
    if (! slice_with_queue[i]) continue;
    slice_ewma_time_[i] = (1-beta_) * slice_ewma_time_[i];
    if (i == slice_id) {
      slice_ewma_time_[i] += beta_ * 1;
    }
  }
  return slice_id;
}

void DownlinkNVSScheduler::SelectFlowsToSchedule (int slice_serve)
{
#ifdef SCHEDULER_DEBUG
	std::cout << "\t Select Flows to schedule" << std::endl;
#endif

  ClearUsersToSchedule();
  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  std::fill(slice_priority_.begin(), slice_priority_.end(), 0);
  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
	{
	  //SELECT FLOWS TO SCHEDULE
	  RadioBearer *bearer = (*it);
    int user_id = bearer->GetUserID();

    if (user_to_slice_[user_id] != slice_serve)
      continue;

	  if (bearer->HasPackets() && bearer->GetDestination ()->GetNodeState () == NetworkNode::STATE_ACTIVE) {
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
		  ENodeB::UserEquipmentRecord *ueRecord = 
        enb->GetUserEquipmentRecord (bearer->GetDestination ()->GetIDNetworkNode ());
		  std::vector<double> spectralEfficiency;
		  std::vector<int> cqiFeedbacks = ueRecord->GetCQI ();
		  int numberOfCqi = cqiFeedbacks.size ();
      AMCModule *amc = GetMacEntity()->GetAmcModule();
		  for (int i = 0; i < numberOfCqi; i++) {
        spectralEfficiency.push_back(amc->GetEfficiencyFromCQI(cqiFeedbacks.at(i)));
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
    if (is_nongreedy_)
      RBsAllocationNonGreedyPF();
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
          << " cumu_bytes: " << user->m_bearers[i]->GetCumulateBytes()
          << " cumu_rbs: " << user->m_bearers[i]->GetCumulateRBs()
          << " hol_delay: " << user->m_bearers[i]->GetHeadOfLinePacketDelay()
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
  for (auto it = users->begin(); it != users->end(); it++) {
    UserToSchedule *ue = *it;
    if (ue->GetListOfAllocatedRBs()->size() > 0) {
      std::vector<double> estimatedSinrValues;

      std::cout << "User(" << ue->GetUserID() << ") allocated RBGS:";
      for (size_t i = 0; i < ue->GetListOfAllocatedRBs()->size (); i++ ) {
        int rbid = ue->GetListOfAllocatedRBs()->at(i);
        if (rbid % rbg_size == 0)
          std::cout << " " << rbid / rbg_size << "(" << ue->GetCqiFeedbacks().at(rbid) << ")";
        
        double sinr = amc->GetSinrFromCQI(
          ue->GetCqiFeedbacks().at(
          ue->GetListOfAllocatedRBs()->at(i)));
        estimatedSinrValues.push_back (sinr);
      }
      double effectiveSinr = GetEesmEffectiveSinr(estimatedSinrValues);
      std::cout << " final_cqi: " << amc->GetCQIFromSinr(effectiveSinr) << std::endl;
      int mcs = amc->GetMCSFromCQI(amc->GetCQIFromSinr(effectiveSinr));
      int transportBlockSize = amc->GetTBSizeFromMCS(mcs, ue->GetListOfAllocatedRBs()->size());

      #if defined(FIRST_SYNTHETIC_EXP) || defined(SECOND_SYNTHETIC_EXP)
      // calculate the transport block size as if the user can have multiple mcs
      transportBlockSize = 0;
      for (int i = 0; i < estimatedSinrValues.size(); i++) {
         transportBlockSize += amc->GetTBSizeFromMCS(amc->GetMCSFromCQI(amc->GetCQIFromSinr(estimatedSinrValues[i])), 1);
      }
      #endif
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
  double averageRate = 1;
  int slice_id = user_to_slice_[user->GetUserID()];
  for (int i = 0; i < MAX_BEARERS; ++i) {
    if (user->m_bearers[i]) {
      averageRate += user->m_bearers[i]->GetAverageTransmissionRate();
    }
  }
  spectralEfficiency = spectralEfficiency * 180000 / 1000; // set the unit to kbps
  averageRate /= 1000.0; // set the unit of both params to kbps
  SchedulerAlgoParam param = slice_algo_params_[slice_id];
  if (param.alpha == 0) {
    metric = pow(spectralEfficiency, param.epsilon) / pow(averageRate, param.psi);
  }
  else {
    // the prioritized flow has no packet, set metric to 0
    if (user->m_dataToTransmit[slice_priority_[slice_id]] == 0) {
      metric = 0;
    }
    else {
      RadioBearer* bearer = user->m_bearers[slice_priority_[slice_id]];
      double HoL = bearer->GetHeadOfLinePacketDelay();
      metric = HoL * pow(spectralEfficiency, param.epsilon)
          / pow(averageRate, param.psi);
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

void
DownlinkNVSScheduler::RBsAllocationNonGreedyPF() {
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
