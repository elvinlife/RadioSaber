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

#include "dl-pf-packet-scheduler.h"
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
#include "../../../core/idealMessages/ideal-control-messages.h"
#include <fstream>

DL_PF_PacketScheduler::DL_PF_PacketScheduler(std::string config_fname="")
{
  std::ifstream ifs(config_fname, std::ifstream::in);
  if (ifs.is_open()) {
    int begin_id = 0;
    ifs >> num_slices_;

    for (int i = 0; i < num_slices_; ++i)
      ifs >> slice_weights_[i];
    for (int i = 0; i < num_slices_; ++i) {
      slice_exp_time_[i] = slice_weights_[i];
      int num_ue;
      ifs >> num_ue;
      for (int j = 0; j < num_ue; ++j) {
        user_to_slice_[begin_id + j] = i;
      }
      begin_id += num_ue;
    }
    for (int i = 0; i < num_slices_; ++i) {
      int algo;
      ifs >> algo;
      if (algo == 0)
        slice_algo_[i] = MT;
      else if (algo == 1)
        slice_algo_[i] = PF;
      else if (algo == 2)
        slice_algo_[i] = MLWDF;
      else
        slice_algo_[i] = PF;
    }
  }
  else {
    throw std::runtime_error("Fail to open configuration file");
  }
  ifs.close();
  SetMacEntity (0);
  CreateFlowsToSchedule ();
}

DL_PF_PacketScheduler::~DL_PF_PacketScheduler()
{
  Destroy ();
}

void
DL_PF_PacketScheduler::DoStopSchedule (void)
{
#ifdef SCHEDULER_DEBUG
  std::cout << "\t Creating Packet Burst" << std::endl;
#endif

  PacketBurst* pb = new PacketBurst ();

  //Create Packet Burst
  FlowsToSchedule *flowsToSchedule = GetFlowsToSchedule ();

  for (FlowsToSchedule::iterator it = flowsToSchedule->begin (); it != flowsToSchedule->end (); it++)
    {
	  FlowToSchedule *flow = (*it);

	  int availableBytes = flow->GetAllocatedBits ()/8;

	  if (availableBytes > 0)
	    {
		  flow->GetBearer()->UpdateTransmittedBytes (availableBytes);
      flow->GetBearer()->UpdateCumulateRBs (flow->GetListOfAllocatedRBs()->size());
      int app_id = flow->GetBearer()->GetApplication()->GetApplicationID();

      std::cerr << GetTimeStamp()
          << " flow: " << app_id
          << " cumu_bytes: " << flow->GetBearer()->GetCumulateBytes()
          << " cumu_rbs: " << flow->GetBearer()->GetCumulateRBs()
          << " hol_delay: " << flow->GetBearer()->GetHeadOfLinePacketDelay()
          << " user: " << app_id
          << " slice: " << user_to_slice_[app_id]
          << std::endl;

	      RlcEntity *rlc = flow->GetBearer ()->GetRlcEntity ();
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
	  else
	    {}
    }
  UpdateTimeStamp();

  //SEND PACKET BURST

#ifdef SCHEDULER_DEBUG
  if (pb->GetNPackets () == 0)
    std::cout << "\t Send only reference symbols" << std::endl;
#endif

  GetMacEntity ()->GetDevice ()->SendPacketBurst (pb);
}

double
DL_PF_PacketScheduler::ComputeSchedulingMetric (RadioBearer *bearer, double spectralEfficiency, int subChannel)
{
  /*
   * For the PF scheduler the metric is computed
   * as follows:
   *
   * metric = spectralEfficiency / averageRate
   */
  double metric = 0;
  // if (intra_sched_ == PF) {
  //   metric = (spectralEfficiency * 180000.)
	//  				  / bearer->GetAverageTransmissionRate();
  // }
  // else if (intra_sched_ == MT) {
  //   metric = spectralEfficiency;
  // }
  // else if (intra_sched_ == MLWDF) {
  //   if (bearer->GetApplication()->GetApplicationType() == 
  //     Application::APPLICATION_TYPE_INFINITE_BUFFER) {
  //       metric = (spectralEfficiency * 180000.)
	//  				  / bearer->GetAverageTransmissionRate();
  //     }
  //   else {
  //     double HOL = bearer->GetHeadOfLinePacketDelay();
  //     metric = HOL * (spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate();
  //   }
  // }
  metric = (spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate();
  return metric;
}

