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

#ifndef DOWNLINKNVSSCHEDULER_H_
#define DOWNLINKNVSSCHEDULER_H_

#include "packet-scheduler.h"

class DownlinkNVSScheduler: public PacketScheduler {
	enum Scheduler {MT, PF, TTA, MLWDF};
private:
	int		user_to_slice_[MAX_APPS];
	double	slice_weights_[MAX_SLICES];
	double 	slice_exp_time_[MAX_SLICES];

	int		num_slices_ = 1;
	int		schedule_scheme_ = 1;

	const double	beta_ = 0.1;
	Scheduler		intra_sched_;

public:
	DownlinkNVSScheduler(std::string config_fname="");
	virtual ~DownlinkNVSScheduler();

	void SelectSliceToServe( int& );
	void SelectFlowsToSchedule ( int );

	virtual void DoSchedule (void);
	virtual void DoStopSchedule (void);

	virtual void RBsAllocation ();
	void RBsAllocationForUE();
	virtual double ComputeSchedulingMetric (
		RadioBearer* bearer, double spectralEfficiency,
		int subChannel );
	virtual double ComputeSchedulingMetric (
		UserToSchedule* user, double spectralEfficiency );

	void UpdateAverageTransmissionRate (int);

};

#endif /* DOWNLINKPACKETSCHEDULER_H_ */
