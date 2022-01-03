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

#ifndef DOWNLINKORACLESCHEDULER_H_
#define DOWNLINKORACLESCHEDULER_H_

#include "packet-scheduler.h"

class DownlinkSubOptScheduler: public PacketScheduler {
	enum Scheduler {MT, PF, TTA};
private:
	int		type1_app_[MAX_APPS];
	int		type1_bitrates_[MAX_SLICES];

	int		type2_app_[MAX_APPS];
	double	type2_weights_[MAX_SLICES];
	double	type2_rbs_offset_[MAX_SLICES];

	int		num_type2_slices_		= 1;
	int		num_type1_slices_		= 1;
	int		num_type1_apps_			= 0;

	const double beta_			= 0.1;
	
	Scheduler intra_sched_ = PF;
	


public:
	DownlinkSubOptScheduler(std::string config_fname);
	virtual ~DownlinkSubOptScheduler();

	void SelectFlowsToSchedule ();

	virtual void DoSchedule (void);
	virtual void DoStopSchedule (void);

	virtual void RBsAllocation ();
	virtual double ComputeSchedulingMetric (RadioBearer *bearer,
			double spectralEfficiency,
			int subChannel,
			double widebandEfficiency);

	void UpdateAverageTransmissionRate (void);
};

#endif /* DOWNLINKPACKETSCHEDULER_H_ */
