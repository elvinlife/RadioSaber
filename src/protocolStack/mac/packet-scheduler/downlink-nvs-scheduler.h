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
	enum Scheduler {MT, PF, TTA};
private:
	int		type1_app_[MAX_APPS];
	int		type2_app_[MAX_APPS];
	int		type1_bitrates_[MAX_SLICES];
	int		type1_exp_bitrates_[MAX_SLICES];
	double	type2_weights_[MAX_SLICES];
	double 	type2_exp_time_[MAX_SLICES];

	int		num_type2_slices_		= 1;
	int		num_type1_slices_		= 1;
	int		num_type1_apps_			= 0;

	const double beta_			= 0.1;

	Scheduler intra_sched_ = PF;



public:
	DownlinkNVSScheduler(std::string config_fname="");
	virtual ~DownlinkNVSScheduler();

	void SelectSliceToServe(int&, bool&);
	void SelectFlowsToSchedule ();

	virtual void DoSchedule (void);
	virtual void DoStopSchedule (void);

	virtual void RBsAllocation ();
	virtual double ComputeSchedulingMetric (
		RadioBearer *bearer, double spectralEfficiency,
		int subChannel, double wideBandEfficiency);

	void UpdateAverageTransmissionRate (void);

	bool isType1(int app_id) {
		return app_id < num_type1_apps_;
	}
};

#endif /* DOWNLINKPACKETSCHEDULER_H_ */
