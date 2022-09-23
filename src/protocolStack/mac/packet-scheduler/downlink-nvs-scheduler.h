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

#ifndef DOWNLINKNVSSCHEDULER_H_
#define DOWNLINKNVSSCHEDULER_H_

#include "packet-scheduler.h"
#include <vector>

class DownlinkNVSScheduler: public PacketScheduler {
	enum Scheduler {MT, PF, TTA, MLWDF};
private:
	int       user_to_slice_[MAX_UES];
	double    slice_weights_[MAX_SLICES];
	double    slice_exp_time_[MAX_SLICES];
	Scheduler slice_algo_[MAX_SLICES];
  
	int       slice_priority_[MAX_SLICES];
	int       num_slices_ = 0;
	int       schedule_scheme_ = 1;
  bool      is_nongreedy_;

  // the beta_ for inter-slice scheduling
  const double beta_ = 0.01;

public:
	DownlinkNVSScheduler(std::string config_fname="", bool is_nongreedy = false);
	virtual ~DownlinkNVSScheduler();

	int  SelectSliceToServe();
	void SelectFlowsToSchedule ( int );

	virtual void DoSchedule (void);
	virtual void DoStopSchedule (void);

	virtual void RBsAllocation ();
	virtual double ComputeSchedulingMetric (
		UserToSchedule* user, double spectralEfficiency );
	void UpdateAverageTransmissionRate (int);

  void RBsAllocationOptimalPF();
  double AssignRBsGivenMCS(std::vector<int>& assigned_mcs,
      std::vector<UserToSchedule*>& rbgs_assignment);
};

#endif /* DOWNLINKPACKETSCHEDULER_H_ */
