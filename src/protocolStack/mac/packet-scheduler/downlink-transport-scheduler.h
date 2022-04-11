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
 */

#ifndef DOWNLINKTRANSPORTSCHEDULER_H_
#define DOWNLINKTRANSPORTSCHEDULER_H_

#include "packet-scheduler.h"

class DownlinkTransportScheduler: public PacketScheduler {
	enum Scheduler {MT, PF, TTA, MLWDF};
private:
	int       user_to_slice_[MAX_APPS];
	double    slice_weights_[MAX_SLICES];
	Scheduler slice_algo_[MAX_SLICES];
	int       slice_priority_[MAX_SLICES];
	double    slice_rbs_offset_[MAX_SLICES];

	int       num_slices_ = 1;
	int       schedule_scheme_ = 1;

	const double beta_ = 0.1;
	int			inter_sched_ = 0;

public:
	DownlinkTransportScheduler(std::string config_fname, int algo);
	virtual ~DownlinkTransportScheduler();

	void SelectFlowsToSchedule ();

	virtual void DoSchedule (void);
	virtual void DoStopSchedule (void);

	virtual void RBsAllocation ();
	void FinalizeAllocation();

	virtual double ComputeSchedulingMetric (
		UserToSchedule* user,
		double spectralEfficiency
	);

	void UpdateAverageTransmissionRate (void);
};

#endif /* DOWNLINKPACKETSCHEDULER_H_ */
