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


#include "downlink-transport-scheduler.h"
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
#include <utility>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <unordered_map>
#include <limits>
using std::vector;
using std::unordered_map;

using coord_t = std::pair<int, int>;
using coord_cqi_t = std::pair<coord_t, double>;

DownlinkTransportScheduler::DownlinkTransportScheduler(std::string config_fname, int algo)
{
  std::ifstream ifs(config_fname);
  if (ifs.is_open()) {
    int begin_id = 0;
    ifs >> schedule_scheme_ >> num_slices_;

    for (int i = 0; i < num_slices_; ++i)
      ifs >> slice_weights_[i];
    for (int i = 0; i < num_slices_; ++i) {
      slice_rbs_offset_[i] = 0;
      int num_ue;
      ifs >> num_ue;
      for (int j = 0; j < num_ue; ++j) {
        user_to_slice_[begin_id + j] = i;
      }
      begin_id += num_ue;
    }
  }
  else {
    throw std::runtime_error("Fail to open configuration file");
  }
  ifs.close();
  inter_sched_ = algo;
  SetMacEntity (0);
  CreateUsersToSchedule();
}

DownlinkTransportScheduler::~DownlinkTransportScheduler()
{
  Destroy ();
}


void DownlinkTransportScheduler::SelectFlowsToSchedule ()
{
#ifdef SCHEDULER_DEBUG
	std::cout << "\t Select Flows to schedule" << std::endl;
#endif

  //ClearFlowsToSchedule ();
  ClearUsersToSchedule();

  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  //std::cerr << GetTimeStamp();
  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
	{
	  //SELECT FLOWS TO SCHEDULE
	  RadioBearer *bearer = (*it);

	  if (bearer->HasPackets() && bearer->GetDestination ()->GetNodeState () == NetworkNode::STATE_ACTIVE)
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
      AMCModule *amc = GetMacEntity()->GetAmcModule();
		  for (int i = 0; i < numberOfCqi; i++)
			{
        spectralEfficiency.push_back(amc->GetEfficiencyFromCQI(cqiFeedbacks.at(i)));
			}
		  //create flow to scheduler record

      InsertFlowToUser(bearer, dataToTransmit, spectralEfficiency, cqiFeedbacks);
		}
	  else
	    {}
	}
  //std::cerr << std::endl;
}

void
DownlinkTransportScheduler::DoSchedule (void)
{
#ifdef SCHEDULER_DEBUG
	std::cout << "Start DL packet scheduler for node "
			<< GetMacEntity ()->GetDevice ()->GetIDNetworkNode()<< std::endl;
#endif

  UpdateAverageTransmissionRate ();
  SelectFlowsToSchedule ();

  if (GetUsersToSchedule()->size() == 0)
	{}
  else
	{
	  //RBsAllocation ();
    RBsAllocationForUE();
	}

  StopSchedule ();
}

void
DownlinkTransportScheduler::DoStopSchedule(void)
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
              << " user: " << user->GetUserID()
              << " flow: " << user->m_bearers[i]->GetApplication()->GetApplicationID()
              << " cumu_bytes: " << user->m_bearers[i]->GetCumulateBytes()
              << " cumu_rbs: " << user->m_bearers[i]->GetCumulateRBs()
              << " hol_delay: " << user->m_bearers[i]->GetHeadOfLinePacketDelay()
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

static unordered_map<int, vector<int>> UpperBound(double** flow_spectraleff, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  double sum_bits = 0;
  unordered_map<int, vector<int>> slice_rbgs;
  for (int j = 0; j < nb_slices; ++j) {
    if (slice_quota_rbgs[j] <= 0)
      continue;
    vector<std::pair<int, double>> sorted_cqi;
    for (int i = 0; i < nb_rbgs; ++i) {
      sorted_cqi.emplace_back(i, flow_spectraleff[i][j]);
    }
    std::sort(sorted_cqi.begin(), sorted_cqi.end(),
      [](std::pair<int, double> a, std::pair<int, double> b) {
        return a.second > b.second;
        }
      );
    for (int k = 0; k < slice_quota_rbgs[j]; ++k) {
      slice_rbgs[j].push_back(sorted_cqi[k].first);
      sum_bits += sorted_cqi[k].second;
    }
  }
  fprintf(stderr, "all_bytes: %.0f\n", sum_bits * 180 / 8 * 4); // 4 for rbg_size
  return slice_rbgs;
}

// return rb id to slice id mapping vector
static vector<int> GreedyByRow(double** flow_spectraleff, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  vector<int> slice_rbgs(nb_slices, 0);
  vector<int> rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i) {
    double max_eff = -1;
    int assigned_slice = -1;
    for (int j = 0; j < nb_slices; ++j) {
      if (flow_spectraleff[i][j] > max_eff && slice_rbgs[j] < slice_quota_rbgs[j]) {
        max_eff = flow_spectraleff[i][j];
        assigned_slice = j;
      }
    }
    assert(assigned_slice != -1);
    rbg_to_slice[i] = assigned_slice;
    slice_rbgs[assigned_slice] += 1;
  }
  double sum_bits = 0;
  for (int i = 0; i < nb_rbgs; ++i) {
    sum_bits += flow_spectraleff[i][rbg_to_slice[i]];
  }
  fprintf(stderr, "all_bytes: %.0f\n", sum_bits * 180 / 8 * 4); // 4 for rbg_size
  return rbg_to_slice;
}

static vector<int> SubOpt(double** flow_spectraleff, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  vector<int> slice_rbgs(nb_slices, 0);
  vector<int> rbg_to_slice(nb_rbgs, -1);
  // set negative to 0 first
  for (int i = 0; i < slice_quota_rbgs.size(); ++i) {
    if (slice_quota_rbgs[i] < 0)
      slice_quota_rbgs[i] = 0;
  }
  for (int i = 0; i < nb_rbgs; ++i) {
    double max_eff = -1;
    int assigned_slice = -1;
    for (int j = 0; j < nb_slices; ++j) {
      if (flow_spectraleff[i][j] > max_eff) {
        max_eff = flow_spectraleff[i][j];
        assigned_slice = j;
      }
    }
    assert(assigned_slice != -1);
    rbg_to_slice[i] = assigned_slice;
    slice_rbgs[assigned_slice] += 1;
  }
  // do the reallocation
  unordered_map<int, int> slice_more;
  unordered_map<int, int> slice_fewer;
  for(int i = 0; i < nb_slices; ++i) {
    if (slice_rbgs[i] > slice_quota_rbgs[i]) {
      slice_more[i] = slice_rbgs[i] - slice_quota_rbgs[i];
    }
    else if (slice_rbgs[i] < slice_quota_rbgs[i]) {
      slice_fewer[i] = slice_quota_rbgs[i] - slice_rbgs[i];
    }
  }
  while (slice_more.size() > 0 && slice_fewer.size() > 0) {
    // one reallocation
    //std::cout << slice_more.size() << "*" << slice_fewer.size() << std::endl;
    int from_slice = -1, to_slice, rbg_id;
    double min_tbs_reduce = std::numeric_limits<double>::max();
    for (int i = 0; i < nb_rbgs; ++i) {
      if ( slice_more.find(rbg_to_slice[i]) != slice_more.end() ) {
        for (auto it = slice_fewer.begin(); it != slice_fewer.end(); ++it) {
          assert(flow_spectraleff[i][rbg_to_slice[i]] >= flow_spectraleff[i][it->first]);
          if (flow_spectraleff[i][rbg_to_slice[i]] - flow_spectraleff[i][it->first] < min_tbs_reduce) {
            min_tbs_reduce = flow_spectraleff[i][rbg_to_slice[i]] - flow_spectraleff[i][it->first];
            from_slice = rbg_to_slice[i];
            to_slice = it->first;
            rbg_id = i;
          }
        }
      }
    }
    assert(from_slice != -1);
    slice_rbgs[from_slice] -= 1;
    slice_rbgs[to_slice] += 1;
    rbg_to_slice[rbg_id] = to_slice;
    slice_more.at(from_slice) -= 1;
    slice_fewer.at(to_slice) -= 1;
        // std::cout << "reallocate: " 
        // << from_slice << " -> " << to_slice << " "
        // << slice_rbgs[from_slice] << " & " << slice_rbgs[to_slice] << " "
        // << slice_more[from_slice] << " & " << slice_fewer[to_slice] << " "
        // << std::endl;
    if (slice_more.at(from_slice) <= 0 || slice_rbgs[from_slice] <= 0) {
      slice_more.erase(from_slice);
    }
    if (slice_fewer.at(to_slice) <= 0) {
      slice_fewer.erase(to_slice);
    }
  }
  double sum_bits = 0;
  for (int i = 0; i < nb_rbgs; ++i) {
    sum_bits += flow_spectraleff[i][rbg_to_slice[i]];
  }
  fprintf(stderr, "all_bytes: %.0f\n", sum_bits * 180 / 8 * 4); // 4 for rbg_size
  return rbg_to_slice;
}

static vector<int> MaximizeCell(double** flow_spectraleff, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  // it's ok that the quota is negative
  vector<coord_cqi_t> sorted_cqi;
  vector<int>         slice_rbgs(nb_slices, 0);
  vector<int>         rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i)
    for (int j = 0; j < nb_slices; ++j) {
      sorted_cqi.emplace_back(coord_t(i, j), flow_spectraleff[i][j]);
    }
    std::sort(sorted_cqi.begin(), sorted_cqi.end(), [](coord_cqi_t a, coord_cqi_t b){return a.second > b.second; });
    for (auto it = sorted_cqi.begin(); it != sorted_cqi.end(); ++it) {
      int rbg_id = it->first.first;
      int slice_id = it->first.second;
      if (slice_rbgs[slice_id] < slice_quota_rbgs[slice_id] && rbg_to_slice[rbg_id] == -1) {
        rbg_to_slice[rbg_id] = slice_id;
        slice_rbgs[slice_id] += 1;
    } 
  }
  double sum_bits = 0;
  for (int i = 0; i < nb_rbgs; ++i) {
    sum_bits += flow_spectraleff[i][rbg_to_slice[i]];
  }
  fprintf(stderr, "all_bytes: %.0f\n", sum_bits * 180 / 8 * 4); // 4 for rbg_size
  return rbg_to_slice;
}

static vector<int> VogelApproximate(double** flow_spectraleff, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  vector<int>   slice_rbgs(nb_slices, 0);
  vector<int>   rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i) {
    int max_diff = -1;
    coord_t coord_1st;
    coord_t coord_2nd;
    // horizonal search
    for (int j = 0; j < nb_rbgs; ++j) {
      // the rbg has been allocated
      if (rbg_to_slice[j] != -1)
        continue;
      double eff_1st = -1, eff_2nd = -1;
      int slice_1st, slice_2nd;

      for (int k = 0; k < nb_slices; ++k) {
        // the slice has reached quota
        if (slice_rbgs[k] >= slice_quota_rbgs[k])
          continue;
        if (eff_1st == -1 || flow_spectraleff[j][k] > eff_1st) {
          slice_1st = k;
          eff_1st = flow_spectraleff[j][k];
          continue;
        }
        if (eff_2nd == -1 || flow_spectraleff[j][k] > eff_2nd) {
          slice_2nd = k;
          eff_2nd = flow_spectraleff[j][k];
          continue;
        }
      }
      if (eff_1st - eff_2nd > max_diff) {
        max_diff = eff_1st - eff_2nd;
        coord_1st = coord_t(j, slice_1st);
        coord_2nd = coord_t(j, slice_2nd);
      }
    }
    // vertical search
    for (int k = 0; k < nb_slices; ++k) {
      // the slice has reached quota
      if (slice_rbgs[k] >= slice_quota_rbgs[k])
        continue;
      double eff_1st = -1, eff_2nd = -1;
      int rbg_1st, rbg_2nd;
      for (int j = 0; j < nb_rbgs; ++j) {
        if (rbg_to_slice[j] != -1)
          continue;
        if (eff_1st == -1 || flow_spectraleff[j][k] > eff_1st) {
          rbg_1st = j;
          eff_1st = flow_spectraleff[j][k];
          continue;
        }
        if (eff_2nd == -1 || flow_spectraleff[j][k] > eff_2nd) {
          rbg_2nd = j;
          eff_2nd = flow_spectraleff[j][k];
          continue;
        }
      }
      if (eff_1st - eff_2nd > max_diff) {
        max_diff = eff_1st - eff_2nd;
        coord_1st = coord_t(rbg_1st, k);
        coord_2nd = coord_t(rbg_2nd, k);
      }
    }
    rbg_to_slice[coord_1st.first] = coord_1st.second;
    slice_rbgs[coord_1st.second] += 1;
  }
  double sum_bits = 0;
  for (int i = 0; i < nb_rbgs; ++i) {
    sum_bits += flow_spectraleff[i][rbg_to_slice[i]];
  }
  fprintf(stderr, "all_bytes: %.0f\n", sum_bits * 180 / 8 * 4); // 4 for rbg_size
  return rbg_to_slice;
}

void
DownlinkTransportScheduler::RBsAllocationForUE()
{
  UsersToSchedule* users = GetUsersToSchedule();
  int nb_rbs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  int rbg_size = get_rbg_size(nb_rbs);
  // currently nb_rbgs should be divisible
  nb_rbs = nb_rbs - (nb_rbs % rbg_size);
  assert(nb_rbs % rbg_size == 0);

  // find out slices without data/flows at all, and assign correct rb target
  std::vector<bool> slice_with_data(users->size(), false);
  std::vector<int> slice_target_rbs(users->size(), 0);
  int num_nonempty_slices = 0;
  int extra_rbs = nb_rbs;
  for (auto it = users->begin(); it != users->end(); ++it) {
    int user_id = (*it)->GetUserID();
    int slice_id = user_to_slice_[user_id];
    if (slice_with_data[slice_id])
      continue;
    num_nonempty_slices += 1;
    slice_with_data[slice_id] = true;
    slice_target_rbs[slice_id] = (int)(nb_rbs * slice_weights_[slice_id] + slice_rbs_offset_[slice_id]);
    extra_rbs -= slice_target_rbs[slice_id];
  }
  assert(num_nonempty_slices != 0);
  // we enable reallocation between slices, but not flows
  bool is_first_slice = true;
  int rand_begin_idx = rand();
  for (int i = 0; i < num_slices_; ++i) {
    int k = (i + rand_begin_idx) % num_slices_;
    if (slice_with_data[k]) {
      slice_target_rbs[k] += extra_rbs / num_nonempty_slices;
      if (is_first_slice) {
        slice_target_rbs[k] += extra_rbs % num_nonempty_slices;
        is_first_slice = false;
      }
    }
  }
  int nb_rbgs = nb_rbs / rbg_size;
  // calculate the rbg quota for slices
  std::vector<int> slice_quota_rbgs(num_slices_, 0);
  std::vector<int> slice_final_rbgs(num_slices_, 0);
  int extra_rbgs = nb_rbgs;
  for (int i = 0; i < num_slices_; ++i) {
    slice_quota_rbgs[i] = (int)(slice_target_rbs[i] / rbg_size);
    extra_rbgs -= slice_quota_rbgs[i];
  }
  is_first_slice = true;
  rand_begin_idx = rand();
  for (int i = 0; i < num_slices_; ++i) {
    int k = (rand_begin_idx + i) % num_slices_;
    if(slice_with_data[k]) {
      slice_quota_rbgs[k] += extra_rbgs / num_nonempty_slices;
      if (is_first_slice) {
        slice_quota_rbgs[k] += extra_rbgs % num_nonempty_slices;
        is_first_slice = false;
      }
    }
  }

  std::cout << "slice target RBs:";
  for (int i = 0; i < num_slices_; ++i) {
    std::cout << "(" << i << ", " << slice_target_rbs[i] << ", " << slice_quota_rbgs[i] << ") ";
  }
  std::cout << std::endl;

  // create a matrix of flow metrics (RBG, flow index)
  double metrics[nb_rbgs][users->size ()];

  for (int i = 0; i < nb_rbgs; i++) {
    for (int j = 0; j < users->size(); j++) {
      metrics[i][j] = ComputeSchedulingMetric(
        users->at(j),
        users->at(j)->GetSpectralEfficiency().at(i * rbg_size)
      );
    }
  }

  AMCModule *amc = GetMacEntity ()->GetAmcModule ();

  // 1st index: rbg_id; 2nd index: slice_id
  // the flow_id when the rbg is assigned to the slice
  int **user_index = new int*[nb_rbgs];
  double **flow_spectraleff = new double*[nb_rbgs];

  // Assign the flow_id and cqi to every slice in every rbg
  for (int i = 0; i < nb_rbgs; i++)
  {
    user_index[i] = new int[num_slices_];
    flow_spectraleff[i] = new double[num_slices_];
    for (int k = 0; k < num_slices_; k++) {
      user_index[i][k] = -1;
      flow_spectraleff[i][k] = 0;
    }
    std::vector<double> max_ranks(num_slices_, -1);     // the highest flow metric in every slice
    for (int j = 0; j < users->size(); ++j) {
      int user_id = users->at(j)->GetUserID();
      int slice_id = user_to_slice_[user_id];
      if (metrics[i][j] > max_ranks[slice_id]) {
        max_ranks[slice_id] = metrics[i][j];
        user_index[i][slice_id] = j;
        flow_spectraleff[i][slice_id] = users->at(j)->GetSpectralEfficiency().at(i * rbg_size);
      }
    }
  }

  // calculate the assignment of rbgs to slices
  vector<int> rbg_to_slice;
  unordered_map<int, vector<int>> slice_rbgs;
  if (inter_sched_ == 0) {
    rbg_to_slice = GreedyByRow(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else if (inter_sched_ == 1) {
    rbg_to_slice = SubOpt(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else if (inter_sched_ == 2) {
    rbg_to_slice = MaximizeCell(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else if (inter_sched_ == 3 ) {
    rbg_to_slice = VogelApproximate(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else {
    slice_rbgs = UpperBound(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }

  // ToDo: Generalize the framework
  if (inter_sched_ < 4 ) {
    for (int i = 0; i < rbg_to_slice.size(); ++i) {
      int uindex = user_index[i][rbg_to_slice[i]];
      assert(uindex != -1);
      int sid = user_to_slice_[users->at(uindex)->GetUserID()];
      assert(sid == rbg_to_slice[i]);
      slice_final_rbgs[sid] += 1;
      int l = i * rbg_size, r = (i+1) * rbg_size;
      for (int j = l; j < r; ++j) {
        users->at(uindex)->GetListOfAllocatedRBs()->push_back(j);
      }
    }
  }
  else {
    for (auto it = slice_rbgs.begin(); it != slice_rbgs.end(); ++it) {
      int sid = it->first;
      auto rbg_list = it->second;
      slice_final_rbgs[sid] += rbg_list.size();
      for (int i = 0; i < rbg_list.size(); ++i) {
        int uindex = user_index[rbg_list[i]][sid];
        assert(uindex != -1);
        int l = rbg_list[i] * rbg_size, r = (rbg_list[i] + 1) * rbg_size;
        for (int j = l; j < r; ++j) {
          users->at(uindex)->GetListOfAllocatedRBs()->push_back(j);
        }
      }
    }
  }

  for (int i = 0; i < num_slices_; ++i) {
    slice_rbs_offset_[i] = slice_target_rbs[i] - slice_final_rbgs[i] * rbg_size;
  }

  // free the flow_id and flow_spectraleff
  for (int i = 0; i < nb_rbgs; i++) {
    delete[] user_index[i];
    delete[] flow_spectraleff[i];
  }
  delete[] user_index;
  delete[] flow_spectraleff;

  PdcchMapIdealControlMessage *pdcchMsg = new PdcchMapIdealControlMessage();
  for (auto it = users->begin(); it != users->end(); it++) {
    UserToSchedule *ue = *it;
    if (ue->GetListOfAllocatedRBs()->size() > 0) {
      std::vector<double> estimatedSinrValues;
      for (int rb = 0; rb < ue->GetListOfAllocatedRBs()->size (); rb++ ) {
        double sinr = amc->GetSinrFromCQI (
        ue->GetCqiFeedbacks ().at (
          ue->GetListOfAllocatedRBs ()->at (rb)));
        estimatedSinrValues.push_back (sinr);
      }
      double effectiveSinr = GetEesmEffectiveSinr(estimatedSinrValues);
      int mcs = amc->GetMCSFromCQI(amc->GetCQIFromSinr(effectiveSinr));
      int transportBlockSize = amc->GetTBSizeFromMCS(mcs, ue->GetListOfAllocatedRBs()->size());
      ue->UpdateAllocatedBits(transportBlockSize);
      for (int rb = 0; rb < ue->GetListOfAllocatedRBs()->size(); rb++) {
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

void
DownlinkTransportScheduler::RBsAllocation ()
{
#ifdef SCHEDULER_DEBUG
	std::cout << " ---- DownlinkTransportScheduler::RBsAllocation";
#endif
  FlowsToSchedule* flows = GetFlowsToSchedule ();
  int nb_rbs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  int rbg_size = get_rbg_size(nb_rbs);
  // currently nb_rbgs should be divisible
  nb_rbs = nb_rbs - (nb_rbs % rbg_size);
  assert(nb_rbs % rbg_size == 0);

  // find out slices without data/flows at all, and assign correct rb target
  std::vector<bool> slice_with_data(num_slices_, false);
  std::vector<int> slice_target_rbs(num_slices_, 0);
  int num_nonempty_slices = 0;
  int extra_rbs = nb_rbs;
  for (int k = 0; k < flows->size(); k++)
  {
    int app_id = flows->at(k)->GetBearer()->GetApplication()->GetApplicationID();
    int slice_id = user_to_slice_[app_id];
    if (slice_with_data[slice_id])
      continue;
    num_nonempty_slices += 1;
    slice_with_data[slice_id] = true;
    slice_target_rbs[slice_id] = (int)(nb_rbs * slice_weights_[slice_id] + slice_rbs_offset_[slice_id]);
    extra_rbs -= slice_target_rbs[slice_id];
  }
  assert(num_nonempty_slices != 0);
  // we enable reallocation between slices, but not flows
  bool is_first_slice = true;
  int rand_begin_idx = rand();
  for (int i = 0; i < num_slices_; ++i) {
    int k = (i + rand_begin_idx) % num_slices_;
    if (slice_with_data[k]) {
      slice_target_rbs[k] += extra_rbs / num_nonempty_slices;
      if (is_first_slice) {
        slice_target_rbs[k] += extra_rbs % num_nonempty_slices;
        is_first_slice = false;
      }
    }
  }
  int nb_rbgs = nb_rbs / rbg_size;
  // calculate the rbg quota for slices
  std::vector<int> slice_quota_rbgs(num_slices_, 0);
  std::vector<int> slice_final_rbgs(num_slices_, 0);
  int extra_rbgs = nb_rbgs;
  for (int i = 0; i < num_slices_; ++i) {
    slice_quota_rbgs[i] = (int)(slice_target_rbs[i] / rbg_size);
    extra_rbgs -= slice_quota_rbgs[i];
  }
  is_first_slice = true;
  rand_begin_idx = rand();
  for (int i = 0; i < num_slices_; ++i) {
    int k = (rand_begin_idx + i) % num_slices_;
    if(slice_with_data[k]) {
      slice_quota_rbgs[k] += extra_rbgs / num_nonempty_slices;
      if (is_first_slice) {
        slice_quota_rbgs[k] += extra_rbgs % num_nonempty_slices;
        is_first_slice = false;
      }
    }
  }

  std::cout << "slice target RBs:";
  for (int i = 0; i < num_slices_; ++i) {
    std::cout << "(" << i << ", " << slice_target_rbs[i] << ", " << slice_quota_rbgs[i] << ") ";
  }
  std::cout << std::endl;

  // create a matrix of flow metrics (RBG, flow index)
  double metrics[nb_rbgs][flows->size ()];

  for (int i = 0; i < nb_rbgs; i++) {
	  for (int j = 0; j < flows->size (); j++) {
		  metrics[i][j] = ComputeSchedulingMetric (
        flows->at (j)->GetBearer (),
        flows->at (j)->GetSpectralEfficiency ().at (i * rbg_size),
        i);
	  }
  }

#ifdef SCHEDULER_DEBUG
  //std::cout << ", available RBGs " << nb_rbgs << ", flows " << flows->size () << std::endl;
  for (int ii = 0; ii < flows->size (); ii++)
  {
	std::cout << "\t metrics for flow "
		  << flows->at (ii)->GetBearer ()->GetApplication ()->GetApplicationID () << ":";
	for (int jj = 0; jj < nb_rbgs; jj++) {
      fprintf(stdout, " (%d, %.3f, %d, %.3f)",
          jj, metrics[jj][ii], 
          flows->at(ii)->GetCqiFeedbacks().at(jj * rbg_size),
          flows->at(ii)->GetSpectralEfficiency().at(jj * rbg_size));
	  }
	std::cout << std::endl;
  }
#endif

  AMCModule *amc = GetMacEntity ()->GetAmcModule ();

  // 1st index: rbg_id; 2nd index: slice_id
  // the flow_id when the rbg is assigned to the slice
  int **flow_id = new int*[nb_rbgs];
  double **flow_spectraleff = new double*[nb_rbgs];

  // Assign the flow_id and cqi to every slice in every rbg
  for (int s = 0; s < nb_rbgs; s++)
  {
    flow_id[s] = new int[num_slices_];
    flow_spectraleff[s] = new double[num_slices_];
    for (int k = 0; k < num_slices_; k++) {
      flow_id[s][k] = -1;
      flow_spectraleff[s][k] = 0;
    }

    std::vector<double> max_ranks(num_slices_, -1);     // the highest flow metric in every slice
    for (int k = 0; k < flows->size(); k++)
    {
      int app_id = flows->at(k)->GetBearer()->GetApplication()->GetApplicationID();
      int slice_id = user_to_slice_[app_id];
      if (metrics[s][k] > max_ranks[slice_id]) {
        max_ranks[slice_id] = metrics[s][k];
        flow_id[s][slice_id] = k;
        flow_spectraleff[s][slice_id] = flows->at(k)->GetSpectralEfficiency().at(s * rbg_size);
      }
    }
  }

  // print the slice_quota and slice_cqi for debugging
  // if (GetTimeStamp() < 1000) {
  //   fprintf(stderr, "slice quota: ");
  //   for (int i = 0; i < num_type2_slices_; i++) {
  //     fprintf(stderr, "%d ", slice_quota_rbgs[i]);
  //   }
  //   fprintf(stderr, "\n");
  //   for (int i = 0; i < nb_rbgs; i++) {
  //     for (int j = 0; j < num_type2_slices_; j++) {
  //       fprintf(stderr, "(%d, %.0f) ", flow_id[i][j], flow_spectraleff[i][j] * 180 / 8 * 4);
  //     }
  //     fprintf(stderr, "\n");
  //   }
  // }

  // calculate the assignment of rbgs to slices
  vector<int> rbg_to_slice;
  unordered_map<int, vector<int>> slice_rbgs;
  if (inter_sched_ == 0) {
    rbg_to_slice = GreedyByRow(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else if (inter_sched_ == 1) {
    rbg_to_slice = SubOpt(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else if (inter_sched_ == 2) {
    rbg_to_slice = MaximizeCell(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else if (inter_sched_ == 3 ) {
    rbg_to_slice = VogelApproximate(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }
  else {
    slice_rbgs = UpperBound(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
  }

  // ToDo: Generalize the framework
  if (inter_sched_ < 4 ) {
    for (int i = 0; i < rbg_to_slice.size(); ++i) {
      int fid = flow_id[i][rbg_to_slice[i]];
      assert(fid != -1);
      int sid = user_to_slice_[flows->at(fid)->GetBearer()->GetApplication()->GetApplicationID()];
      assert(sid == rbg_to_slice[i]);
      slice_final_rbgs[sid] += 1;
      int l = i * rbg_size, r = (i+1) * rbg_size;
      for (int j = l; j < r; ++j) {
        flows->at(fid)->GetListOfAllocatedRBs()->push_back(j);
      }
    }
  }
  else {
    for (auto it = slice_rbgs.begin(); it != slice_rbgs.end(); ++it) {
      int sid = it->first;
      auto rbg_list = it->second;
      slice_final_rbgs[sid] += rbg_list.size();
      for (int i = 0; i < rbg_list.size(); ++i) {
        int fid = flow_id[rbg_list[i]][sid];
        assert(fid != -1);
        int l = rbg_list[i] * rbg_size, r = (rbg_list[i] + 1) * rbg_size;
        for (int j = l; j < r; ++j) {
          flows->at(fid)->GetListOfAllocatedRBs()->push_back(j);
        }
      }
    }
  }

  for (int i = 0; i < num_slices_; ++i) {
    slice_rbs_offset_[i] = slice_target_rbs[i] - slice_final_rbgs[i] * rbg_size;
  }

  // free the flow_id and flow_spectraleff
  for (int i = 0; i < nb_rbgs; i++) {
    delete[] flow_id[i];
    delete[] flow_spectraleff[i];
  }
  delete[] flow_id;
  delete[] flow_spectraleff;

  FinalizeAllocation();  
}

void
DownlinkTransportScheduler::FinalizeAllocation()
{
  FlowsToSchedule* flows = GetFlowsToSchedule ();
  AMCModule *amc = GetMacEntity ()->GetAmcModule ();
  int nb_rbs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  int rbg_size = get_rbg_size(nb_rbs);

  //Finalize the allocation
  PdcchMapIdealControlMessage *pdcchMsg = new PdcchMapIdealControlMessage ();

  for (FlowsToSchedule::iterator it = flows->begin (); it != flows->end (); it++)
  {
    FlowToSchedule *flow = (*it);

    if (flow->GetListOfAllocatedRBs ()->size () > 0)
    {
	    std::cout << "Flow(" << flow->GetBearer()->GetApplication()->GetApplicationID() << ")";
	    for (int rb = 0; rb < flow->GetListOfAllocatedRBs()->size(); rb++) {
        int rbid = flow->GetListOfAllocatedRBs()->at(rb);
        if (rbid % rbg_size == 0)
          std::cout << " " << rbid / rbg_size;
	    }
	    std::cout << std::endl;

      //this flow has been scheduled
      std::vector<double> estimatedSinrValues;
      for (int rb = 0; rb < flow->GetListOfAllocatedRBs ()->size (); rb++ ) {
          double sinr = amc->GetSinrFromCQI (
                  flow->GetCqiFeedbacks ().at (flow->GetListOfAllocatedRBs ()->at (rb)));
          estimatedSinrValues.push_back (sinr);
      }
      //compute the effective sinr
      double effectiveSinr = GetEesmEffectiveSinr (estimatedSinrValues);
      //get the MCS for transmission
      int mcs = amc->GetMCSFromCQI (amc->GetCQIFromSinr (effectiveSinr));
      //define the amount of bytes to transmit
      int transportBlockSize = amc->GetTBSizeFromMCS (mcs, flow->GetListOfAllocatedRBs ()->size ());

      // double effectiveSinr = estimatedSinrValues[0];
      // int mcs = amc->GetMCSFromCQI(amc->GetCQIFromSinr(effectiveSinr));
      // int transportBlockSize = 0;
      // for (int i = 0; i < estimatedSinrValues.size(); i++) {
      //   transportBlockSize += amc->GetTBSizeFromMCS(
      //     amc->GetMCSFromCQI(
      //       amc->GetCQIFromSinr(
      //         estimatedSinrValues[i])), 1);
      // }

      flow->UpdateAllocatedBits (transportBlockSize);
#ifdef SCHEDULER_DEBUG
		  std::cout << "\t\t --> flow "	<< flow->GetBearer ()->GetApplication ()->GetApplicationID ()
				  << " has been scheduled: " <<
				  "\n\t\t\t nb of RBs " << flow->GetListOfAllocatedRBs ()->size () <<
				  "\n\t\t\t effectiveSinr " << effectiveSinr <<
				  "\n\t\t\t tbs " << transportBlockSize <<
          "\n\t\t\t data to transmit " << flow->GetDataToTransmit() <<
				  "\n\t\t\t cqi " << amc->GetCQIFromSinr(effectiveSinr)
				  << std::endl;
#endif
		//create PDCCH messages
		for (int rb = 0; rb < flow->GetListOfAllocatedRBs ()->size (); rb++ )
		{
		  pdcchMsg->AddNewRecord (PdcchMapIdealControlMessage::DOWNLINK,
				  flow->GetListOfAllocatedRBs ()->at (rb),
								  flow->GetBearer ()->GetDestination (),
								  mcs);
		}
	  }
  }

  if (pdcchMsg->GetMessage()->size () > 0)
  {
    GetMacEntity ()->GetDevice ()->GetPhy ()->SendIdealControlMessage (pdcchMsg);
  }
  delete pdcchMsg;
}

double
DownlinkTransportScheduler::ComputeSchedulingMetric(RadioBearer *bearer, double spectralEfficiency, int subChannel)
{
  double metric = 0;
  switch (intra_sched_) {
    case MT:
      metric = spectralEfficiency;
      break;
    case PF:
      metric = (spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate();
      break;
    case MLWDF:
    {
      double a = -log10 (0.05) / 1;
      double HOL = bearer->GetHeadOfLinePacketDelay ();
      metric = (a * HOL) * ((spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate ());
      break;
    }
    default:
      metric = spectralEfficiency;
  }
  return metric;
}

double
DownlinkTransportScheduler::ComputeSchedulingMetric(UserToSchedule* user, double spectralEfficiency)
{
  double metric = 0;
  double averageRate = 1;
  for (int i = 0; i < MAX_BEARERS; ++i) {
    if (user->m_bearers[i]) {
      averageRate += user->m_bearers[i]->GetAverageTransmissionRate();
    }
  }
  if (schedule_scheme_ == 0) {
    switch (intra_sched_) {
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
    if (user->m_dataToTransmit[GetHighestPriority()] == 0){
      metric = 0;
    }
    else {
      if (GetHighestPriority() == 1) {
        RadioBearer* bearer = user->m_bearers[1];
        double HoL = bearer->GetHeadOfLinePacketDelay();
        metric = HoL * (spectralEfficiency * 180000.) / averageRate;
      }
      else if (GetHighestPriority() == 0) {
        metric = (spectralEfficiency * 180000.) / averageRate;
      }
      else {
        throw std::runtime_error(
          "Invalid bearer priority: " + std::to_string(GetHighestPriority()));
      }
    }
  }
  return metric;
}

void
DownlinkTransportScheduler::UpdateAverageTransmissionRate (void)
{
  // we should update the user average transmission rate instead of the flow transmission rate
  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
    {
      RadioBearer *bearer = (*it);
      bearer->UpdateAverageTransmissionRate ();
    }
}