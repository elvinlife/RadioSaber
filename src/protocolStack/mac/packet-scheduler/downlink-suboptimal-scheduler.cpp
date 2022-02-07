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
 *         Lukasz Rajewski <lukasz.rajewski@gmail.com> (optimized PRB allocation)
 */


#include "downlink-suboptimal-scheduler.h"
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
#include "../../../flows/MacQueue.h"
#include "../../../utility/eesm-effective-sinr.h"
#include <cstdio>
#include <utility>
#include <unordered_map>
#include <fstream>
#include <cassert>

DownlinkSubOptScheduler::DownlinkSubOptScheduler(std::string config_fname)
{
  std::ifstream ifs(config_fname);
  if (ifs.is_open()) {
    int begin_id = 0;
    ifs >> num_type1_slices_ >> num_type2_slices_;
    for (int i = 0; i < num_type1_slices_; ++i)
      ifs >> type1_bitrates_[i];
    for (int i = 0; i < num_type1_slices_; ++i) {
      int num_ue;
      ifs >> num_ue;
      for (int j = 0; j < num_ue; ++j) {
        type1_app_[begin_id + j] = i;
      }
      begin_id += num_ue;
    }
    num_type1_apps_ = begin_id;

    for (int i = 0; i < num_type2_slices_; ++i)
      ifs >> type2_weights_[i];
    for (int i = 0; i < num_type2_slices_; ++i) {
      type2_rbs_offset_[i] = 0;
      int num_ue;
      ifs >> num_ue;
      for (int j = 0; j < num_ue; ++j) {
        type2_app_[begin_id + j] = i;
      }
      begin_id += num_ue;
    }
  }
  else {
    throw std::runtime_error("Fail to open configuration file");
  }
  ifs.close();
  
  SetMacEntity (0);
  CreateFlowsToSchedule ();
}

DownlinkSubOptScheduler::~DownlinkSubOptScheduler()
{
  Destroy ();
}

void DownlinkSubOptScheduler::SelectFlowsToSchedule ()
{
#ifdef SCHEDULER_DEBUG
	std::cout << "\t Select Flows to schedule" << std::endl;
#endif

  ClearFlowsToSchedule ();

  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  //std::cerr << GetTimeStamp();
  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
	{
	  //SELECT FLOWS TO SCHEDULE
	  RadioBearer *bearer = (*it);

	  if (bearer->HasPackets () && bearer->GetDestination ()->GetNodeState () == NetworkNode::STATE_ACTIVE)
		{
		  //compute data to transmit
		  int dataToTransmit;
		  if (bearer->GetApplication ()->GetApplicationType () == Application::APPLICATION_TYPE_INFINITE_BUFFER)
			{
			  dataToTransmit = 100000000;
			}
		  else
			{
			  dataToTransmit = bearer->GetQueueSize ();
			}

		  //compute spectral efficiency
		  ENodeB *enb = (ENodeB*) GetMacEntity ()->GetDevice ();
		  ENodeB::UserEquipmentRecord *ueRecord = enb->GetUserEquipmentRecord (bearer->GetDestination ()->GetIDNetworkNode ());
		  std::vector<double> spectralEfficiency;
		  std::vector<int> cqiFeedbacks = ueRecord->GetCQI ();
		  int numberOfCqi = cqiFeedbacks.size ();
      AMCModule *amc = GetMacEntity()->GetAmcModule();
		  for (int i = 0; i < numberOfCqi; i++)
			{
        spectralEfficiency.push_back(amc->GetEfficiencyFromCQI(cqiFeedbacks.at(i)));
			}
		  //create flow to scheduler record
      if (dataToTransmit > 0)
		    InsertFlowToSchedule(bearer, dataToTransmit, spectralEfficiency, cqiFeedbacks);
		}
	  else
	    {}
	}
  //std::cerr << std::endl;
}

void
DownlinkSubOptScheduler::DoSchedule (void)
{
#ifdef SCHEDULER_DEBUG
	std::cout << "Start DL packet scheduler for node "
			<< GetMacEntity ()->GetDevice ()->GetIDNetworkNode()<< std::endl;
#endif

  UpdateAverageTransmissionRate ();
  SelectFlowsToSchedule ();

  if (GetFlowsToSchedule ()->size() == 0)
	{}
  else
	{
	  RBsAllocation ();
	}

  StopSchedule ();
}

void
DownlinkSubOptScheduler::DoStopSchedule (void)
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

		  flow->GetBearer()->UpdateTransmittedBytes (min(availableBytes, flow->GetDataToTransmit()));
      flow->GetBearer()->UpdateCumulateRBs (flow->GetListOfAllocatedRBs()->size());

      std::cerr << GetTimeStamp()
          << " flow: " << flow->GetBearer()->GetApplication()->GetApplicationID()
          << " cumu_bytes: " << flow->GetBearer()->GetCumulateBytes()
          << " cumu_rbs: " << flow->GetBearer()->GetCumulateRBs()
          << " hol_delay: " << flow->GetBearer()->GetHeadOfLinePacketDelay()
          << std::endl;

	      RlcEntity *rlc = flow->GetBearer ()->GetRlcEntity ();
	      PacketBurst* pb2 = rlc->TransmissionProcedure (availableBytes);

// #ifdef SCHEDULER_DEBUG
// 	      std::cout << "\t\t  nb of packets: " << pb2->GetNPackets () << std::endl;
// #endif

	      if (pb2->GetNPackets () > 0)
	        {
	    	  std::list<Packet*> packets = pb2->GetPackets ();
	    	  std::list<Packet* >::iterator it;
	    	  for (it = packets.begin (); it != packets.end (); it++)
	    	    {
// #ifdef SCHEDULER_DEBUG
// 	    		  std::cout << "\t\t  added packet of bytes " << (*it)->GetSize () << std::endl;
// #endif
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
  //UpdateAverageTransmissionRate ();

  //SEND PACKET BURST

#ifdef SCHEDULER_DEBUG
  if (pb->GetNPackets () == 0)
    std::cout << "\t Send only reference symbols" << std::endl;
#endif

  GetMacEntity ()->GetDevice ()->SendPacketBurst (pb);
}

void
DownlinkSubOptScheduler::RBsAllocation ()
{
#ifdef SCHEDULER_DEBUG
	std::cout << " ---- DownlinkSubOptScheduler::RBsAllocation";
#endif
  FlowsToSchedule* flows = GetFlowsToSchedule ();
  int nbOfRBs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  int rbg_size = get_rbg_size(nbOfRBs);
  int nbOfGroups = (nbOfRBs + rbg_size - 1) / rbg_size;

  std::vector<double> sliceTargetRBs(num_type2_slices_);
  std::cout << "\t slice quota:";
  for (int i = 0; i < num_type2_slices_; ++i) {
    double quota = nbOfRBs * type2_weights_[i];
    sliceTargetRBs[i] = quota + type2_rbs_offset_[i];
    std::cout << "\t" << sliceTargetRBs[i];
  }
  std::cout << std::endl;

  // create a matrix of flow metrics (RBG, flow index)
  double metrics[nbOfGroups][flows->size ()];
  // create a matrix of slices information
  int sliceFlow[nbOfGroups][num_type2_slices_];
  double sliceMaxMetric[nbOfGroups][num_type2_slices_];
  double sliceCQI[nbOfGroups][num_type2_slices_];
  // initialize
  for (int i = 0; i < nbOfGroups; ++i) {
    for (int j = 0; j < num_type2_slices_; ++j) {
      sliceFlow[i][j] = sliceMaxMetric[i][j] = sliceCQI[i][j] = -1;
    }
  }

  int rbgToSlice[nbOfGroups];
  std::vector<int> sliceRBs(num_type2_slices_, 0);

  for (int i = 0; i < nbOfGroups; i++) {
	  for (int j = 0; j < flows->size (); j++) {
      int app_id = flows->at(j)->GetBearer()->GetApplication()->GetApplicationID();
      int slice_id = type2_app_[app_id];
		  metrics[i][j] = ComputeSchedulingMetric (
        flows->at (j)->GetBearer (),
        flows->at (j)->GetSpectralEfficiency ().at (i * rbg_size),
        i, flows->at(j)->GetAllEfficiency());
      if (metrics[i][j] > sliceMaxMetric[i][slice_id]) {
        sliceMaxMetric[i][slice_id] = metrics[i][j];
        sliceFlow[i][slice_id] = j;
        sliceCQI[i][slice_id] = flows->at(j)->GetSpectralEfficiency().at(i * rbg_size);
      }
	  }
    
  }
  for (int i = 0; i < nbOfGroups; ++i) {
    double maxCQI = 0;
    for (int j = 0; j < num_type2_slices_; ++j) {
      //fprintf(stderr, "rb: %d slice: %d cqi: %.3f\n", i, j, sliceCQI[i][j]);
      if (sliceCQI[i][j] > maxCQI) {
        rbgToSlice[i] = j;
        maxCQI = sliceCQI[i][j];
      }
    }
    sliceRBs[rbgToSlice[i]] += rbg_size;
  }

  // the reallocation procedure
  // fprintf(stderr, "reallocation begin\n");
  std::unordered_map<int, double> slice_more;
  std::unordered_map<int, double> slice_fewer;
  for (int i = 0; i < num_type2_slices_; ++i) {
    if (sliceRBs[i] >= sliceTargetRBs[i] + rbg_size) {
      slice_more[i] = sliceRBs[i] - (sliceTargetRBs[i] > 0 ? sliceTargetRBs[i] : 0);
    }
    else if (sliceRBs[i] <= sliceTargetRBs[i] - rbg_size) {
      slice_fewer[i] = sliceTargetRBs[i] - sliceRBs[i];
    }
  }
  while (slice_more.size() > 0 && slice_fewer.size() > 0) {
    std::cout << "more: ";
    for (auto it = slice_more.begin(); it != slice_more.end(); ++it) {
        std::cout << it->first << ", ";
    }
    std::cout << std::endl;
    for (auto it = slice_fewer.begin(); it != slice_fewer.end(); ++it) {
        std::cout << it->first << ", ";
    }
    std::cout << std::endl;
    double min_tbs_reduce = 100;
    int rb_id, from_slice = -1, to_slice;
    for (int i = 0; i < nbOfGroups; ++i) {
      if ( slice_more.find(rbgToSlice[i]) != slice_more.end() ) {
        for (auto it = slice_fewer.begin(); it != slice_fewer.end(); ++it) {
          if (sliceCQI[i][rbgToSlice[i]] < sliceCQI[i][it->first]) {
            char error_log[200];
            sprintf(error_log, "rb: %d optimal_s: %d %d abnormal_s: %.3f %.3f",
                i, rbgToSlice[i], it->first, 
                sliceCQI[i][rbgToSlice[i]], sliceCQI[i][it->first]);
            throw std::runtime_error(std::string(error_log));
          } 
          if (sliceCQI[i][rbgToSlice[i]] - sliceCQI[i][it->first] < min_tbs_reduce) {
            min_tbs_reduce = sliceCQI[i][rbgToSlice[i]] - sliceCQI[i][it->first];
            rb_id = i;
            from_slice = rbgToSlice[i];
            to_slice = it->first;
          }
        }
      }
    }
    if (from_slice == -1)
      break;
    sliceRBs[from_slice] -= rbg_size;
    sliceRBs[to_slice] += rbg_size;
    rbgToSlice[rb_id] = to_slice;
    slice_more.at(from_slice) -= rbg_size;
    slice_fewer.at(to_slice) -= rbg_size;
    std::cout << "reallocate: " 
        << from_slice << " -> " << to_slice << " "
        << sliceRBs[from_slice] << " & " << sliceRBs[to_slice] << " "
        << slice_more[from_slice] << " & " << slice_fewer[to_slice] << " "
        << "rbg_size: " << rbg_size
        << std::endl;
    if (slice_more.at(from_slice) <= 0 || sliceRBs[from_slice] <= 0) {
      slice_more.erase(from_slice);
    }
    if (slice_fewer.at(to_slice) <= 0) {
      slice_fewer.erase(to_slice);
    }
  }
  

#ifdef SCHEDULER_DEBUG
  // for (int ii = 0; ii < flows->size (); ii++)
  //   {
	//   std::cout << "\t metrics for flow "
	// 		  << flows->at (ii)->GetBearer ()->GetApplication ()->GetApplicationID () << ":";
	//   for (int jj = 0; jj < nbOfGroups; jj++)
	//     {
  //       fprintf(stdout, " (%d, %.3f, %d)",
  //           jj, metrics[jj][ii], 
  //           flows->at(ii)->GetCqiFeedbacks().at(jj * rbg_size));
	//     }
	//   std::cout << std::endl;
  //   }
#endif

  AMCModule *amc = GetMacEntity ()->GetAmcModule ();
  bool * l_bFlowScheduled = new bool[flows->size ()];
  int l_iScheduledFlows = 0;
  std::vector<double> * l_bFlowScheduledSINR = new std::vector<double>[flows->size ()];
  for (int k = 0; k < flows->size (); k++)
      l_bFlowScheduled[k] = false;

  //RBs allocation
  for (int s = 0; s < nbOfGroups; s++)
    {
      bool SubbandAllocated = true;
      FlowToSchedule* scheduledFlow;
      int sliceServe = rbgToSlice[s];
      int l_iScheduledFlowIndex = sliceFlow[s][sliceServe];
      scheduledFlow = flows->at(l_iScheduledFlowIndex);

      if (SubbandAllocated)
        {
          // allocate the sth rbg
          int l = s*rbg_size, r = (s+1)*rbg_size;
          if (r > nbOfRBs) r = nbOfRBs;

          for (int i = l; i < r; ++i) {
            scheduledFlow->GetListOfAllocatedRBs()->push_back(i);
            double sinr = amc->GetSinrFromCQI(scheduledFlow->GetCqiFeedbacks().at(i));
            l_bFlowScheduledSINR[l_iScheduledFlowIndex].push_back(sinr);
          }

          double effectiveSinr = GetEesmEffectiveSinr (l_bFlowScheduledSINR[l_iScheduledFlowIndex]);
          int mcs = amc->GetMCSFromCQI (amc->GetCQIFromSinr (effectiveSinr));
          int alloc_num = scheduledFlow->GetListOfAllocatedRBs()->size();
          int transportBlockSize = amc->GetTBSizeFromMCS (mcs, alloc_num);
          if (transportBlockSize >= scheduledFlow->GetDataToTransmit() * 8)
          {
              l_bFlowScheduled[l_iScheduledFlowIndex] = true;
              l_iScheduledFlows++;
          }
        }
    }

  for (int i = 0; i < num_type2_slices_; ++i) {
    type2_rbs_offset_[i] = sliceTargetRBs[i] - sliceRBs[i];
  }

  delete [] l_bFlowScheduled;
  delete [] l_bFlowScheduledSINR;

  //Finalize the allocation
  PdcchMapIdealControlMessage *pdcchMsg = new PdcchMapIdealControlMessage ();

  for (FlowsToSchedule::iterator it = flows->begin (); it != flows->end (); it++)
    {
      FlowToSchedule *flow = (*it);

	  std::cout << "Flow(" << flow->GetBearer()->GetApplication()->GetApplicationID() << ")";
	  for (int rb = 0; rb < flow->GetListOfAllocatedRBs()->size(); rb++) {
      int rbid = flow->GetListOfAllocatedRBs()->at(rb);
      if (rbid % rbg_size == 0)
        std::cout << " " << rbid / rbg_size;
	  }
	  std::cout << std::endl;

      if (flow->GetListOfAllocatedRBs ()->size () > 0)
        {
          //this flow has been scheduled
          std::vector<double> estimatedSinrValues;
          for (int rb = 0; rb < flow->GetListOfAllocatedRBs ()->size (); rb++ ) {
              double sinr = amc->GetSinrFromCQI (
                      flow->GetCqiFeedbacks ().at (flow->GetListOfAllocatedRBs ()->at (rb)));
              estimatedSinrValues.push_back (sinr);
            }

          //compute the effective sinr
          double effectiveSinr = GetEesmEffectiveSinr (estimatedSinrValues);

          //get the MCS for transmission
          int mcs = amc->GetMCSFromCQI (amc->GetCQIFromSinr (effectiveSinr));

          //define the amount of bytes to transmit
          int transportBlockSize = amc->GetTBSizeFromMCS (mcs, flow->GetListOfAllocatedRBs ()->size ());
          flow->UpdateAllocatedBits (transportBlockSize);

#ifdef SCHEDULER_DEBUG
		  std::cout << "\t\t --> flow "	<< flow->GetBearer ()->GetApplication ()->GetApplicationID ()
				  << " has been scheduled: " <<
				  "\n\t\t\t nb of RBs " << flow->GetListOfAllocatedRBs ()->size () <<
				  "\n\t\t\t effectiveSinr " << effectiveSinr <<
				  "\n\t\t\t tbs " << transportBlockSize <<
          "\n\t\t\t data to transmit " << flow->GetDataToTransmit() <<
				  "\n\t\t\t cqi " << amc->GetCQIFromSinr(effectiveSinr)
				  << std::endl;
#endif

		  //create PDCCH messages
		  for (int rb = 0; rb < flow->GetListOfAllocatedRBs ()->size (); rb++ )
		    {
			  pdcchMsg->AddNewRecord (PdcchMapIdealControlMessage::DOWNLINK,
					  flow->GetListOfAllocatedRBs ()->at (rb),
									  flow->GetBearer ()->GetDestination (),
									  mcs);
		    }
	    }
    }

  if (pdcchMsg->GetMessage()->size () > 0)
    {
      GetMacEntity ()->GetDevice ()->GetPhy ()->SendIdealControlMessage (pdcchMsg);
    }
  delete pdcchMsg;
}

double
DownlinkSubOptScheduler::ComputeSchedulingMetric(RadioBearer *bearer, double spectralEfficiency, int subChannel, double wbEff)
{
  double metric = 0;
  switch (intra_sched_) {
    case MT:
      metric = spectralEfficiency;
      break;
    case PF:
      metric = (spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate();
      break;
    case TTA:
      metric = spectralEfficiency / wbEff;
      break;
    case MLWDF:
    {
      double a = -log10 (0.05) / 1;
      double HOL = bearer->GetHeadOfLinePacketDelay ();
      metric = (a * HOL) * ((spectralEfficiency * 180000.) / bearer->GetAverageTransmissionRate ());
      break;
    }

    default:
      metric = spectralEfficiency;
  }
  return metric;
}

void
DownlinkSubOptScheduler::UpdateAverageTransmissionRate (void)
{
  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
    {
      RadioBearer *bearer = (*it);
      bearer->UpdateAverageTransmissionRate ();
    }
}
