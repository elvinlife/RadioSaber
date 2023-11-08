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

#include "downlink-transport-scheduler.h"
#include <jsoncpp/json/json.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <utility>
#include "../../../core/spectrum/bandwidth-manager.h"
#include "../../../device/ENodeB.h"
#include "../../../device/NetworkNode.h"
#include "../../../flows/MacQueue.h"
#include "../../../flows/application/Application.h"
#include "../../../flows/radio-bearer.h"
#include "../../../load-parameters.h"
#include "../../../phy/lte-phy.h"
#include "../../../protocolStack/mac/AMCModule.h"
#include "../../../protocolStack/rrc/rrc-entity.h"
#include "../../../utility/eesm-effective-sinr.h"
#include "../../packet/Packet.h"
#include "../../packet/packet-burst.h"
#include "../mac-entity.h"

using std::unordered_map;
using std::vector;

// Peter: The coordinates are likely to represent the RBG (i) and the slice index (j).
using coord_t = std::pair<int, int>;
using coord_cqi_t = std::pair<coord_t, double>;

/* d_{u,i}, the instantaneous data rate of UE u for RBG i, which is the
   transmission rate (in bps) per HZ. */
double maxThroughputMetric(DownlinkTransportScheduler::UserToSchedule* user,
                           int index) {
  // user: u, index: i * rbg_size
  return user->GetSpectralEfficiency().at(index);
}

/* d_{u,i} / R_{u}, where R_{u} captures the historical RBG allocation for UE u
   as an exponential weighted moving average of the user's throughput, based on
   its data rate for the RBGs it has been assigned so far. */
double proportionalFairnessMetric(
    DownlinkTransportScheduler::UserToSchedule* user, int index) {
  // user: u, index: i * rbg_size
  double averageRate = 1;
  for (int i = 0; i < MAX_BEARERS; ++i) {
    if (user->m_bearers[i]) {
      averageRate += user->m_bearers[i]->GetAverageTransmissionRate();
    }
  }
  return maxThroughputMetric(user, index) / averageRate;
}

/* D_{u,p} * d_{u,i} / R_{u}, where D_{u,p} is the queuing delay experience by
   packet at the head of the queue corresponding to UE u and priority p. */

// Peter: MLWDF is different from previous enterprise scheduling algorithm,
// because it selects the highest priority first, and if there are multiple UEs
// with the same priority, it then runs this metric.
double mLWDFMetric(DownlinkTransportScheduler::UserToSchedule* user,
                   int index) {
  // user: u, index: i * rbg_size

  // @peter
  // calculate the averateRate across all priority channels
  // get the score of the highest priority channel

  double averageRate = 1;
  int max_priority = -1;
  int selected_bearer = -1;
  for (int i = 0; i < MAX_BEARERS; ++i) {
    if (user->m_bearers[i]) {
      averageRate += user->m_bearers[i]->GetAverageTransmissionRate();
      int priority = user->m_bearers[i]->GetPriority();
      if (priority > max_priority) {
        max_priority = priority;
        selected_bearer = i;
      }
    }
  }

  RadioBearer* bearer = user->m_bearers[selected_bearer];
  double HoL = bearer->GetHeadOfLinePacketDelay();
  return HoL * maxThroughputMetric(user, index) / averageRate;
}

// Peter: Save the score to the log 
void DownlinkTransportScheduler::logScore() {
  std::string filename;
  if (inter_metric_ = &maxThroughputMetric){
    filename = "MaxThroughputMetric.txt";
  } else if (inter_metric_ == &proportionalFairnessMetric){
    filename = "ProportionalFairnessMetric.txt";
  } else if (inter_metric_ == &mLWDFMetric){
    filename = "mLWDFMetric.txt";
  } else {
    throw std::runtime_error("Error: invalid inter-slice metric (objective)");
  }

  std::ofstream logFile(filename, std::ios_base::app);
  if (!logFile.is_open()){
    throw std::runtime_error("Failed to open the score log file");
  }

  for (int i = 0; i < slice_score_.size(); i++) {
    logFile << "Slice Index: " << i << ", Score: " << slice_score_[i] << std::endl;
  }

  logFile.close();
}

// peter: reading in the slice cionfiguration
DownlinkTransportScheduler::DownlinkTransportScheduler(
    std::string config_fname, int interslice_algo, int interslice_metric = 0) {
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
    // peter: get the number of users per slice
    num_ue = ues_per_slice[i].asInt();
    for (int j = 0; j < num_ue; j++) {
      // peter: map user to slice vector, important to differentiate between the i and j here.
      user_to_slice_.push_back(i);
    }
  }
  const Json::Value& slice_schemes = obj["slices"];
  for (int i = 0; i < slice_schemes.size(); i++) {
    int n_slices = slice_schemes[i]["n_slices"].asInt();
    for (int j = 0; j < n_slices; j++) {
      // Peter: this is a useful place to get the enterprise scheduler algorithms parameters.
      slice_weights_.push_back(slice_schemes[i]["weight"].asDouble());
      slice_algo_params_.emplace_back(slice_schemes[i]["algo_alpha"].asInt(),
                                      slice_schemes[i]["algo_beta"].asInt(),
                                      slice_schemes[i]["algo_epsilon"].asInt(),
                                      slice_schemes[i]["algo_psi"].asInt());
    }
  }
  // [peter] for each slice, calculate the priority
  slice_priority_.resize(num_slices_);
  std::fill(slice_priority_.begin(), slice_priority_.end(), 0);
  // [Peter] used later to ensure that the allocations of rbgs overtime is overall fair and accords to the SLA
  slice_rbs_offset_.resize(num_slices_);
  std::fill(slice_rbs_offset_.begin(), slice_rbs_offset_.end(), 0);

  inter_sched_ = interslice_algo;

  // Charlie: set the inter-slice metric (objective)
  switch (interslice_metric) {
    case 0:  // sched 9
      inter_metric_ = &maxThroughputMetric;
      break;
    case 1:  // sched 91
      inter_metric_ = &proportionalFairnessMetric;
      break;
    case 2:  // sched 92
      inter_metric_ = &mLWDFMetric;
      break;
    default:
      throw std::runtime_error("Error: invalid inter-slice metric (objective)");
      break;
  }

  SetMacEntity(0);
  CreateUsersToSchedule();
}

DownlinkTransportScheduler::~DownlinkTransportScheduler() {
  Destroy();
}

// [Peter] For each radio bearer in each UE, get the respective spectral efficiency
// [Peter] update the slice's priority to that of the highest UE's priority in the slice
void DownlinkTransportScheduler::SelectFlowsToSchedule() {
#ifdef SCHEDULER_DEBUG
  std::cout << "\t Select Flows to schedule" << std::endl;
#endif

  ClearUsersToSchedule();
  // Peter: RRC stands for radio resource control, it is the place where the slice context is stored
  RrcEntity* rrc =
      GetMacEntity()->GetDevice()->GetProtocolStack()->GetRrcEntity();
  // Peter: this is unexpected, because according to the paper, radio bearer is managed by individual UE,
  //  but here it seems to suggest that
  // radio bearers are managed by thr RRC
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer();

  // [Peter] Fill the slice priority vector with 0s
  std::fill(slice_priority_.begin(), slice_priority_.end(), 0);
  for (std::vector<RadioBearer*>::iterator it = bearers->begin();
       it != bearers->end(); it++) {
    //SELECT FLOWS TO SCHEDULE
    RadioBearer* bearer = (*it);

    // get the amount of data that the UE wants to transmit
    if (bearer->HasPackets() &&
        bearer->GetDestination()->GetNodeState() == NetworkNode::STATE_ACTIVE) {
      //compute data to transmit
      int dataToTransmit;
      if (bearer->GetApplication()->GetApplicationType() ==
          Application::APPLICATION_TYPE_INFINITE_BUFFER) {
        dataToTransmit = 100000000;
      } else {
        dataToTransmit = bearer->GetQueueSize();
      }

      //compute spectral efficiency
      ENodeB* enb = (ENodeB*)GetMacEntity()->GetDevice();
      // [Peter]Get the previous record of each ue
      ENodeB::UserEquipmentRecord* ueRecord = enb->GetUserEquipmentRecord(
          bearer->GetDestination()->GetIDNetworkNode());
      std::vector<double> spectralEfficiency;
      std::vector<int> cqiFeedbacks = ueRecord->GetCQI();
      // Peter: this is important, each UE has multiple CQI, each for a different frequency
      int numberOfCqi = cqiFeedbacks.size();
      // peter: get the Adaptive Modulation And Coding Scheme, which contains the method to convert CQI into spectral efficiency
      AMCModule* amc = GetMacEntity()->GetAmcModule();
      // [Peter]calculate spectral efficiency based on CQI
      for (int i = 0; i < numberOfCqi; i++) {
        spectralEfficiency.push_back(
            amc->GetEfficiencyFromCQI(cqiFeedbacks.at(i)));
      }

      // create flow to scheduler record
      // [Peter]update the slice priority to be the biggest of all UEs in that slice
      int slice_id = user_to_slice_[bearer->GetUserID()];
      if (bearer->GetPriority() > slice_priority_[slice_id]) {
        slice_priority_[slice_id] = bearer->GetPriority();
      }
      InsertFlowToUser(bearer, dataToTransmit, spectralEfficiency,
                       cqiFeedbacks);
    }
  }
}

// [Peter]update the CQI, call the actual scheduling function
// [Peter]entrance of scheduling
void DownlinkTransportScheduler::DoSchedule(void) {
#ifdef SCHEDULER_DEBUG
  std::cout << "Start DL packet scheduler for node "
            << GetMacEntity()->GetDevice()->GetIDNetworkNode() << std::endl;
#endif
  // [Peter] update each bearer's average transmission rate
  UpdateAverageTransmissionRate();
  // [Peter]update each UE's CQI
  SelectFlowsToSchedule();

  if (GetUsersToSchedule()->size() != 0) {
    RBsAllocation();
  }

  StopSchedule();
}

// peter: actually send out the packets.
void DownlinkTransportScheduler::DoStopSchedule(void) {
  PacketBurst* pb = new PacketBurst();
  UsersToSchedule* uesToSchedule = GetUsersToSchedule();
  for (auto it = uesToSchedule->begin(); it != uesToSchedule->end(); it++) {
    UserToSchedule* user = *it;
    int availableBytes = user->GetAllocatedBits() / 8;
    // let's not reallocate RBs between flows firstly
    for (int i = MAX_BEARERS - 1; i >= 0; i--) {
      if (availableBytes <= 0)
        break;
      if (user->m_dataToTransmit[i] > 0) {
        assert(user->m_bearers[i] != NULL);
        int dataTransmitted = min(availableBytes, user->m_dataToTransmit[i]);
        availableBytes -= dataTransmitted;
        user->m_bearers[i]->UpdateTransmittedBytes(dataTransmitted);
        user->m_bearers[i]->UpdateCumulateRBs(
            user->GetListOfAllocatedRBs()->size());
        std::cerr << GetTimeStamp() << " app: "
                  << user->m_bearers[i]->GetApplication()->GetApplicationID()
                  << " cumu_bytes: " << user->m_bearers[i]->GetCumulateBytes()
                  << " cumu_rbs: " << user->m_bearers[i]->GetCumulateRBs()
                  << " hol_delay: "
                  << user->m_bearers[i]->GetHeadOfLinePacketDelay()
                  << " user: " << user->GetUserID()
                  << " slice: " << user_to_slice_[user->GetUserID()]
                  << std::endl;

        RlcEntity* rlc = user->m_bearers[i]->GetRlcEntity();
        PacketBurst* pb2 = rlc->TransmissionProcedure(dataTransmitted);

        if (pb2->GetNPackets() > 0) {
          std::list<Packet*> packets = pb2->GetPackets();
          std::list<Packet*>::iterator it;
          for (it = packets.begin(); it != packets.end(); it++) {
            Packet* p = (*it);
            pb->AddPacket(p->Copy());
          }
        }
        delete pb2;
      }
    }
  }
  UpdateTimeStamp();

  GetMacEntity()->GetDevice()->SendPacketBurst(pb);
}

// peter: the UpperBound function, as it is written, does not take into account that an RBG can only be allocated to a single slice or user.
static unordered_map<int, vector<int>> UpperBound(double** flow_spectraleff,
                                                  vector<int>& slice_quota_rbgs,
                                                  int nb_rbgs, int nb_slices) {
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
              });
    for (int k = 0; k < slice_quota_rbgs[j]; ++k) {
      slice_rbgs[j].push_back(sorted_cqi[k].first);
      sum_bits += sorted_cqi[k].second;
    }
  }
  fprintf(stderr, "all_bytes: %.0f\n",
          sum_bits * 180 / 8 * 4);  // 4 for rbg_size
  return slice_rbgs;
}

// return rb id to slice id mapping vector
static vector<int> GreedyByRow(double** flow_spectraleff,
                               vector<int>& slice_quota_rbgs, int nb_rbgs,
                               int nb_slices) {
  vector<int> slice_rbgs(nb_slices, 0);
  vector<int> rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i) {
    double max_eff = -1;
    int assigned_slice = -1;
    for (int j = 0; j < nb_slices; ++j) {
      if (flow_spectraleff[i][j] > max_eff &&
          slice_rbgs[j] < slice_quota_rbgs[j]) {
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
  fprintf(stderr, "all_bytes: %.0f\n",
          sum_bits * 180 / 8 * 4);  // 4 for rbg_size
  return rbg_to_slice;
}

static vector<int> SubOpt(double** flow_spectraleff,
                          vector<int>& slice_quota_rbgs, int nb_rbgs,
                          int nb_slices) {
  vector<int> slice_rbgs(nb_slices, 0);
  vector<int> rbg_to_slice(nb_rbgs, -1);
  // set negative to 0 first
  for (size_t i = 0; i < slice_quota_rbgs.size(); ++i) {
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
  for (int i = 0; i < nb_slices; ++i) {
    if (slice_rbgs[i] > slice_quota_rbgs[i]) {
      slice_more[i] = slice_rbgs[i] - slice_quota_rbgs[i];
    } else if (slice_rbgs[i] < slice_quota_rbgs[i]) {
      slice_fewer[i] = slice_quota_rbgs[i] - slice_rbgs[i];
    }
  }
  while (slice_more.size() > 0 && slice_fewer.size() > 0) {
    // one reallocation
    //std::cout << slice_more.size() << "*" << slice_fewer.size() << std::endl;
    int from_slice = -1, to_slice, rbg_id;
    double min_tbs_reduce = std::numeric_limits<double>::max();
    for (int i = 0; i < nb_rbgs; ++i) {
      if (slice_more.find(rbg_to_slice[i]) != slice_more.end()) {
        for (auto it = slice_fewer.begin(); it != slice_fewer.end(); ++it) {
          assert(flow_spectraleff[i][rbg_to_slice[i]] >=
                 flow_spectraleff[i][it->first]);
          if (flow_spectraleff[i][rbg_to_slice[i]] -
                  flow_spectraleff[i][it->first] <
              min_tbs_reduce) {
            min_tbs_reduce = flow_spectraleff[i][rbg_to_slice[i]] -
                             flow_spectraleff[i][it->first];
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
  fprintf(stderr, "all_bytes: %.0f\n",
          sum_bits * 180 / 8 * 4);  // 4 for rbg_size
  return rbg_to_slice;
}

// peter: I think this is the RadioSaber's approach when every enterprise returns the candidate
static vector<int> MaximizeCell(double** flow_spectraleff,
                                vector<int>& slice_quota_rbgs, int nb_rbgs,
                                int nb_slices) {
  // it's ok that the quota is negative
  vector<coord_cqi_t> sorted_cqi;
  vector<int> slice_rbgs(nb_slices, 0);
  vector<int> rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i)
    for (int j = 0; j < nb_slices; ++j) {
      sorted_cqi.emplace_back(coord_t(i, j), flow_spectraleff[i][j]);
    }
  std::sort(sorted_cqi.begin(), sorted_cqi.end(),
            [](coord_cqi_t a, coord_cqi_t b) { return a.second > b.second; });
  for (auto it = sorted_cqi.begin(); it != sorted_cqi.end(); ++it) {
    int rbg_id = it->first.first;
    int slice_id = it->first.second;
    if (slice_rbgs[slice_id] < slice_quota_rbgs[slice_id] &&
        rbg_to_slice[rbg_id] == -1) {
      rbg_to_slice[rbg_id] = slice_id;
      slice_rbgs[slice_id] += 1;
    }
  }
  double sum_bits = 0;
  for (int i = 0; i < nb_rbgs; ++i) {
    sum_bits += flow_spectraleff[i][rbg_to_slice[i]];
  }
  fprintf(stderr, "all_bytes: %.0f\n",
          sum_bits * 180 / 8 * 4);  // 4 for rbg_size
  return rbg_to_slice;
}

// peter: I don't bother to understand this, nor do I think it important
static vector<int> VogelApproximate(double** flow_spectraleff,
                                    vector<int>& slice_quota_rbgs, int nb_rbgs,
                                    int nb_slices) {
  vector<int> slice_rbgs(nb_slices, 0);
  vector<int> rbg_to_slice(nb_rbgs, -1);
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
  fprintf(stderr, "all_bytes: %.0f\n",
          sum_bits * 180 / 8 * 4);  // 4 for rbg_size
  return rbg_to_slice;
}

void DownlinkTransportScheduler::RBsAllocation() {
  UsersToSchedule* users = GetUsersToSchedule();
  int nb_rbs = GetMacEntity()
                   ->GetDevice()
                   ->GetPhy()
                   ->GetBandwidthManager()
                   ->GetDlSubChannels()
                   .size();
  int rbg_size = get_rbg_size(nb_rbs);
  // currently nb_rbgs should be divisible
  nb_rbs = nb_rbs - (nb_rbs % rbg_size);
  assert(nb_rbs % rbg_size == 0);

  // find out slices without data/flows at all, and assign correct rb target
  std::vector<bool> slice_with_data(num_slices_, false);
  std::vector<int> slice_target_rbs(num_slices_, 0);
  int num_nonempty_slices = 0;
  int extra_rbs = nb_rbs;
  for (auto it = users->begin(); it != users->end(); ++it) {
    int user_id = (*it)->GetUserID();
    int slice_id = user_to_slice_[user_id];
    if (slice_with_data[slice_id])
      continue;
    num_nonempty_slices += 1;
    slice_with_data[slice_id] = true;
    slice_target_rbs[slice_id] =
        (int)(nb_rbs * slice_weights_[slice_id] + slice_rbs_offset_[slice_id]);
    extra_rbs -= slice_target_rbs[slice_id];
  }
  // std::cout << "slice target RBs:";
  // for (int i = 0; i < num_slices_; ++i) {
  //   std::cout << "(" << i << ", "
  //     << slice_target_rbs[i] << ", "
  //     << slice_rbs_offset_[i] << ","
  //     << slice_weights_[i] << ") ";
  // }
  // std::cout << std::endl;

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
  // peter: Distribute any remaining RBs (extra_rbs) evenly across the slices with data.
  // peter: If there are leftover RBs after the even distribution, allocate them to the first slice
  // (after randomizing the start index with rand_begin_idx to distribute leftovers randomly).
  is_first_slice = true;
  rand_begin_idx = rand();
  for (int i = 0; i < num_slices_; ++i) {
    int k = (rand_begin_idx + i) % num_slices_;
    if (slice_with_data[k]) {
      slice_quota_rbgs[k] += extra_rbgs / num_nonempty_slices;
      if (is_first_slice) {
        slice_quota_rbgs[k] += extra_rbgs % num_nonempty_slices;
        is_first_slice = false;
      }
    }
  }

  std::cout << "slice_id, target_rbs, quota_rbgs: ";
  for (int i = 0; i < num_slices_; ++i) {
    std::cout << "(" << i << ", " << slice_target_rbs[i] << ", "
              << slice_quota_rbgs[i] << ") ";
  }
  std::cout << std::endl;

  // create a matrix of flow metrics (RBG, flow index)
  // Peter: the metics variable here is very important. It is the enterprise score for each flow and each rbg
  double metrics[nb_rbgs][users->size()];

  for (int i = 0; i < nb_rbgs; i++) {
    for (size_t j = 0; j < users->size(); j++) {
      metrics[i][j] = ComputeSchedulingMetric(
          users->at(j), users->at(j)->GetSpectralEfficiency().at(i * rbg_size));
    }
  }

  AMCModule* amc = GetMacEntity()->GetAmcModule();

  // 1st index: rbg_id; 2nd index: slice_id
  // the flow_id when the rbg is assigned to the slice
  int** user_index = new int*[nb_rbgs];
  double** flow_spectraleff = new double*[nb_rbgs];

  // Assign the flow_id and cqi to every slice in every rbg
  for (int i = 0; i < nb_rbgs; i++) {
    user_index[i] = new int[num_slices_];
    flow_spectraleff[i] = new double[num_slices_];
    for (int k = 0; k < num_slices_; k++) {
      user_index[i][k] = -1;
      flow_spectraleff[i][k] = 0;
    }
    std::vector<double> max_ranks(
        num_slices_, -1);  // the highest flow metric in every slice

    /* Charlie: in this part, we may have different inter-slice scheduling
       metrics (objectives) based on interslice_metric, reflected on the
       inter_metric_ member (function pointer) */
    for (size_t j = 0; j < users->size(); ++j) {
      int user_id = users->at(j)->GetUserID();
      int slice_id = user_to_slice_[user_id];
      if (metrics[i][j] > max_ranks[slice_id]) {
        max_ranks[slice_id] = metrics[i][j];
        user_index[i][slice_id] = j;
        flow_spectraleff[i][slice_id] =
            inter_metric_(users->at(j), i * rbg_size);
      }
    }
  }

  // calculate the assignment of rbgs to slices
  vector<int> rbg_to_slice;
  unordered_map<int, vector<int>> slice_rbgs;

  switch (inter_sched_) {
    case 0:
      rbg_to_slice =
          GreedyByRow(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
      break;
    case 1:
      rbg_to_slice =
          SubOpt(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
      break;
    case 2:
      rbg_to_slice = MaximizeCell(flow_spectraleff, slice_quota_rbgs, nb_rbgs,
                                  num_slices_);
      break;
    case 3:
      rbg_to_slice = VogelApproximate(flow_spectraleff, slice_quota_rbgs,
                                      nb_rbgs, num_slices_);
      break;
    default:
      slice_rbgs =
          UpperBound(flow_spectraleff, slice_quota_rbgs, nb_rbgs, num_slices_);
      break;
  }

  // ToDo: Generalize the framework
  if (inter_sched_ < 4) {
    for (size_t i = 0; i < rbg_to_slice.size(); ++i) {
      // Peter: Update each slice's running score
      slice_score_[rbg_to_slice[i]] += flow_spectraleff[i][rbg_to_slice[i]];

      int uindex = user_index[i][rbg_to_slice[i]];
      assert(uindex != -1);
      int sid = user_to_slice_[users->at(uindex)->GetUserID()];
      assert(sid == rbg_to_slice[i]);
      slice_final_rbgs[sid] += 1;
      int l = i * rbg_size, r = (i + 1) * rbg_size;
      for (int j = l; j < r; ++j) {
        users->at(uindex)->GetListOfAllocatedRBs()->push_back(j);
      }
    }
  } else {
    for (auto it = slice_rbgs.begin(); it != slice_rbgs.end(); ++it) {
      int sid = it->first;
      auto rbg_list = it->second;
      slice_final_rbgs[sid] += rbg_list.size();
      for (size_t i = 0; i < rbg_list.size(); ++i) {
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

  PdcchMapIdealControlMessage* pdcchMsg = new PdcchMapIdealControlMessage();
  std::cout << GetTimeStamp() << std::endl;
  for (auto it = users->begin(); it != users->end(); it++) {
    UserToSchedule* ue = *it;
    if (ue->GetListOfAllocatedRBs()->size() > 0) {
      std::vector<double> estimatedSinrValues;

      std::cout << "User(" << ue->GetUserID() << ") allocated RBGS:";
      for (size_t i = 0; i < ue->GetListOfAllocatedRBs()->size(); i++) {
        int rbid = ue->GetListOfAllocatedRBs()->at(i);
        if (rbid % rbg_size == 0)
          std::cout << " " << rbid / rbg_size << "("
                    << ue->GetCqiFeedbacks().at(rbid) << ")";

        double sinr = amc->GetSinrFromCQI(
            ue->GetCqiFeedbacks().at(ue->GetListOfAllocatedRBs()->at(i)));
        estimatedSinrValues.push_back(sinr);
      }
      double effectiveSinr = GetEesmEffectiveSinr(estimatedSinrValues);
      std::cout << " final_cqi: " << amc->GetCQIFromSinr(effectiveSinr)
                << std::endl;
      int mcs = amc->GetMCSFromCQI(amc->GetCQIFromSinr(effectiveSinr));
      int transportBlockSize =
          amc->GetTBSizeFromMCS(mcs, ue->GetListOfAllocatedRBs()->size());

#if defined(FIRST_SYNTHETIC_EXP) || defined(SECOND_SYNTHETIC_EXP)
      // calculate the transport block size as if the user can have multiple mcs
      transportBlockSize = 0;
      for (int i = 0; i < estimatedSinrValues.size(); i++) {
        transportBlockSize += amc->GetTBSizeFromMCS(
            amc->GetMCSFromCQI(amc->GetCQIFromSinr(estimatedSinrValues[i])), 1);
      }
#endif
      ue->UpdateAllocatedBits(transportBlockSize);
      for (size_t rb = 0; rb < ue->GetListOfAllocatedRBs()->size(); rb++) {
        pdcchMsg->AddNewRecord(PdcchMapIdealControlMessage::DOWNLINK,
                               ue->GetListOfAllocatedRBs()->at(rb),
                               ue->GetUserNode(), mcs);
      }
    }
  }
  if (pdcchMsg->GetMessage()->size() > 0) {
    GetMacEntity()->GetDevice()->GetPhy()->SendIdealControlMessage(pdcchMsg);
  }
  delete pdcchMsg;
}

double DownlinkTransportScheduler::ComputeSchedulingMetric(
    UserToSchedule* user, double spectralEfficiency) {
  double metric = 0;
  double averageRate = 1;
  int slice_id = user_to_slice_[user->GetUserID()];
  for (int i = 0; i < MAX_BEARERS; ++i) {
    if (user->m_bearers[i]) {
      averageRate += user->m_bearers[i]->GetAverageTransmissionRate();
    }
  }
  spectralEfficiency =
      spectralEfficiency * 180000 / 1000;  // set the unit to kbps
  averageRate /= 1000.0;  // set the unit of both params to kbps
  SchedulerAlgoParam param = slice_algo_params_[slice_id];
  if (param.alpha == 0) {
    metric =
        pow(spectralEfficiency, param.epsilon) / pow(averageRate, param.psi);
  } else {
    // the prioritized flow has no packet, set metric to 0
    if (user->m_dataToTransmit[slice_priority_[slice_id]] == 0) {
      metric = 0;
    } else {
      RadioBearer* bearer = user->m_bearers[slice_priority_[slice_id]];
      if (param.beta) {
        double HoL = bearer->GetHeadOfLinePacketDelay();
        metric = HoL * pow(spectralEfficiency, param.epsilon) /
                 pow(averageRate, param.psi);
      } else {
        metric = pow(spectralEfficiency, param.epsilon) /
                 pow(averageRate, param.psi);
      }
    }
  }
  return metric;
}

void DownlinkTransportScheduler::UpdateAverageTransmissionRate(void) {
  // we should update the user average transmission rate instead of the flow transmission rate
  RrcEntity* rrc =
      GetMacEntity()->GetDevice()->GetProtocolStack()->GetRrcEntity();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer();

  for (std::vector<RadioBearer*>::iterator it = bearers->begin();
       it != bearers->end(); it++) {
    RadioBearer* bearer = (*it);
    bearer->UpdateAverageTransmissionRate();
  }
}
