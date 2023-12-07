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

#ifndef DOWNLINKTRANSPORTSCHEDULER_H_
#define DOWNLINKTRANSPORTSCHEDULER_H_

#include <vector>

#include "packet-scheduler.h"

class DownlinkTransportScheduler : public PacketScheduler {
 private:
  // below use customizable scheduler params
  int num_slices_ = 1;
  std::vector<int> user_to_slice_;
  std::vector<double> slice_weights_;
  std::vector<SchedulerAlgoParam> slice_algo_params_;
  std::vector<int> slice_priority_;
  std::vector<double> slice_rbs_offset_;

  const double beta_ = 0.1;
  int inter_sched_ = 0;

 public:
  DownlinkTransportScheduler(std::string config_fname, int algo);
  virtual ~DownlinkTransportScheduler();

  void SelectFlowsToSchedule();

  virtual void DoSchedule(void);
  virtual void DoStopSchedule(void);

  virtual void RBsAllocation();
  virtual double ComputeSchedulingMetric(UserToSchedule* user,
                                         double spectralEfficiency);
  void UpdateAverageTransmissionRate(void);
};

#endif /* DOWNLINKPACKETSCHEDULER_H_ */
