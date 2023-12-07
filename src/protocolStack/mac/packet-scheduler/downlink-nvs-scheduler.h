/* project: RadioSaber; Mode: C++
 * Copyright (c) 2021, 2022, 2023, 2024 University of Illinois Urbana Champaign
 *
 * This file is part of RadioSaber, which is a project built upon LTE-Sim in
 * 2022
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

#ifndef DOWNLINKNVSSCHEDULER_H_
#define DOWNLINKNVSSCHEDULER_H_

#include <vector>

#include "packet-scheduler.h"

class DownlinkNVSScheduler : public PacketScheduler {
 private:
  // below use customizable scheduler params
  int num_slices_ = 1;
  std::vector<int> user_to_slice_;
  std::vector<double> slice_weights_;
  std::vector<SchedulerAlgoParam> slice_algo_params_;
  std::vector<int> slice_priority_;
  std::vector<double> slice_ewma_time_;

  bool is_nongreedy_;
  // the ewma beta for inter-slice scheduling
  const double beta_ = 0.01;

 public:
  DownlinkNVSScheduler(std::string config_fname = "",
                       bool is_nongreedy = false);
  virtual ~DownlinkNVSScheduler();

  int SelectSliceToServe();
  void SelectFlowsToSchedule(int);

  virtual void DoSchedule(void);
  virtual void DoStopSchedule(void);

  virtual void RBsAllocation();
  virtual double ComputeSchedulingMetric(UserToSchedule* user,
                                         double spectralEfficiency);
  void UpdateAverageTransmissionRate(int);

  void RBsAllocationNonGreedyPF();
  double AssignRBsGivenMCS(std::vector<int>& assigned_mcs,
                           std::vector<UserToSchedule*>& rbgs_assignment);
};

#endif /* DOWNLINKPACKETSCHEDULER_H_ */
