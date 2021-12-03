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

class DownlinkOracleScheduler: public PacketScheduler {
	enum Scheduler {MT, PF, TTA};
private:
	int appid_to_slice_[MAX_APPS] 		= {0, 0, 1, 1, 2, 2, 3, 3};
	double slices_weights_[MAX_SLICES]	= {0.25, 0.25, 0.25, 0.25};
	int num_slices_				= 4;
	const double beta_			= 0.1;
	
	Scheduler intra_sched_ = PF;
	
	std::vector<double> slices_exp_times_;
	std::vector<double> slices_rbs_offset_;

public:
	DownlinkOracleScheduler(std::string config_fname);
	virtual ~DownlinkOracleScheduler();

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
