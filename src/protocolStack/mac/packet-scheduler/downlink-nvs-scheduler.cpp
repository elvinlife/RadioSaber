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


#include "downlink-nvs-scheduler.h"
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
#include <limits>
#include <fstream>

DownlinkNVSScheduler::DownlinkNVSScheduler(std::string config_fname)
{
  std::ifstream ifs(config_fname, std::ifstream::in);
  if (ifs.is_open()) {
    int begin_id = 0;
    ifs >> num_type1_slices_ >> num_type2_slices_;
    for (int i = 0; i < num_type1_slices_; ++i)
      ifs >> type1_bitrates_[i];
    for (int i = 0; i < num_type1_slices_; ++i) {
      type1_exp_bitrates_[i] = type1_bitrates_[i];
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
      type2_exp_time_[i] = type2_weights_[i];
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

DownlinkNVSScheduler::~DownlinkNVSScheduler()
{
  Destroy ();
}

void DownlinkNVSScheduler::SelectSliceToServe(int& slice_id, bool& is_type1)
{
  double max_score = 0;
  // for (int i = 0; i < num_type1_slices_; ++i) {
  //   std::cout << i << ": gbr: " << type1_bitrates_[i] << " average: " << type1_exp_bitrates_[i] << std::endl;
  // }
  // currently no type1 slice
  /*
  for (int i = 0; i < num_type1_slices_; ++i) {
    if (type1_exp_bitrates_[i] == 0) {
      slice_id = i;
      is_type1 = true;
      return;
    }
    else {
      double score = (double)type1_bitrates_[i] / type1_exp_bitrates_[i];
      if (score > max_score) {
        max_score = score;
        slice_id = i;
        is_type1 = true;
      }
    }
  }
  */

#ifdef SCHEDULER_DEBUG
  for(int i = 0; i < num_type2_slices_; ++i) {
    std::cout << i << ": weight: " << type2_weights_[i] << " time: " << type2_exp_time_[i] << std::endl;
  }
#endif

  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();
  std::vector<bool> slice_with_queue(num_type2_slices_, false);
  for (auto it = bearers->begin(); it != bearers->end(); it++) {
    RadioBearer *bearer = (*it);
    int app_id = bearer->GetApplication()->GetApplicationID();
    if (bearer->HasPackets () && bearer->GetDestination ()->GetNodeState () == NetworkNode::STATE_ACTIVE)
		{
		  int dataToTransmit;
		  if (bearer->GetApplication ()->GetApplicationType () == Application::APPLICATION_TYPE_INFINITE_BUFFER) {
			  dataToTransmit = 100000000;
			}
		  else {
			  dataToTransmit = bearer->GetQueueSize ();
			}
      if (dataToTransmit > 0) {
        slice_with_queue[type2_app_[app_id]] = true;
      }
    }
  }
  for (int i = 0; i < num_type2_slices_; ++i) {
    if (! slice_with_queue[i]) continue;
    if (type2_exp_time_[i] == 0) {
      slice_id = i;
      is_type1 = false;
      break;
    }
    else {
      double score = type2_weights_[i] / type2_exp_time_[i];
      if (score > max_score) {
        max_score = score;
        slice_id = i;
        is_type1 = false;
      }
    }
  }
  for (int i = 0; i < num_type2_slices_; ++i) {
    if (! slice_with_queue[i]) continue;
    type2_exp_time_[i] = (1-beta_) * type2_exp_time_[i];
    if (i == slice_id) {
      type2_exp_time_[i] += beta_ * 1;
    }
  }
}

void DownlinkNVSScheduler::SelectFlowsToSchedule (int slice_serve)
{
  ClearFlowsToSchedule ();

  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
	{
	  //SELECT FLOWS TO SCHEDULE
	  RadioBearer *bearer = (*it);
    int app_id = bearer->GetApplication()->GetApplicationID();

    // if ( !(isType1(app_id) && is_type1 && type1_app_[app_id] == slice_serve)
    //   && !(!isType1(app_id) && !is_type1 && type2_app_[app_id] == slice_serve) )
    //   continue;
    if (type2_app_[app_id] != slice_serve)
      continue;

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
		  for (int i = 0; i < numberOfCqi; i++)
			{
			  double sEff = GetMacEntity ()->GetAmcModule ()->GetEfficiencyFromCQI (cqiFeedbacks.at (i));
			  spectralEfficiency.push_back (sEff);
			}

		  //create flow to scheduler record
		  InsertFlowToSchedule(bearer, dataToTransmit, spectralEfficiency, cqiFeedbacks);
		}
	  else
	    {}
	}
}

void
DownlinkNVSScheduler::DoSchedule (void)
{
#ifdef SCHEDULER_DEBUG
	std::cout << "\nStart DL packet scheduler for node "
			<< GetMacEntity ()->GetDevice ()->GetIDNetworkNode()
      << " ts: " << GetTimeStamp() << std::endl;
#endif

  // some logic to determine the slice to serve
  int slice_serve = 0;
  bool is_type1 = true;
  SelectSliceToServe(slice_serve, is_type1);
  UpdateAverageTransmissionRate (slice_serve);
  SelectFlowsToSchedule (slice_serve);

  if (GetFlowsToSchedule ()->size() == 0)
	{}
  else
	{
	  RBsAllocation ();
	}

  StopSchedule ();
}

void
DownlinkNVSScheduler::DoStopSchedule (void)
{
  PacketBurst* pb = new PacketBurst ();

  //Create Packet Burst
  FlowsToSchedule *flowsToSchedule = GetFlowsToSchedule ();
  
  bool is_type1 = true;
  int slice_id = -1;
  int aggregate_bits = 0;

  for (FlowsToSchedule::iterator it = flowsToSchedule->begin (); it != flowsToSchedule->end (); it++)
    {
	  FlowToSchedule *flow = (*it);

	  int availableBytes = flow->GetAllocatedBits ()/8;

    if (slice_id == -1) {
      int app_id = flow->GetBearer()->GetApplication()->GetApplicationID();
      if (isType1( app_id ) ){
        slice_id = type1_app_[app_id];
      }
      else {
        slice_id = type2_app_[app_id];
        is_type1 = false;
      }
    }
    aggregate_bits += flow->GetAllocatedBits();

	  if (availableBytes > 0)
	    {
		  flow->GetBearer()->UpdateTransmittedBytes (min(availableBytes, flow->GetDataToTransmit()));
      
      flow->GetBearer()->UpdateCumulateRBs (flow->GetListOfAllocatedRBs()->size());

      std::cerr << GetTimeStamp()
          << " flow: " << flow->GetBearer()->GetApplication()->GetApplicationID()
          << " cumu_bytes: " << flow->GetBearer()->GetCumulateBytes()
          << " cumu_rbs: " << flow->GetBearer()->GetCumulateRBs()
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
DownlinkNVSScheduler::RBsAllocation ()
{
#ifdef SCHEDULER_DEBUG
	std::cout << " ---- DownlinkNVSScheduler::RBsAllocation";
#endif

  FlowsToSchedule* flows = GetFlowsToSchedule ();
  int nbOfRBs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  int rbg_size = get_rbg_size(nbOfRBs);
  int nbOfGroups = (nbOfRBs + rbg_size - 1) / rbg_size;

  // create a matrix of flow metrics
  double metrics[nbOfGroups][flows->size ()];
  for (int i = 0; i < nbOfGroups; i++) {
	  for (int j = 0; j < flows->size (); j++) {
		  metrics[i][j] = ComputeSchedulingMetric (
        flows->at (j)->GetBearer (),
        flows->at (j)->GetSpectralEfficiency ().at (i * rbg_size),
        i,
        flows->at (j)->GetAllEfficiency());
	  }
  }

#ifdef SCHEDULER_DEBUG
  for (int ii = 0; ii < flows->size (); ii++)
    {
	  std::cout << "\t metrics for flow "
			  << flows->at (ii)->GetBearer ()->GetApplication ()->GetApplicationID () << ":";
	  for (int jj = 0; jj < nbOfGroups; jj++)
	  {
      fprintf(stdout, " (%d, %.3f, %d)",
          jj, metrics[jj][ii], 
          flows->at(ii)->GetCqiFeedbacks().at(jj * rbg_size));
	  }
	  std::cout << std::endl;
    }
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
      if (l_iScheduledFlows == flows->size ())
          break;

      double targetMetric = std::numeric_limits<double>::min();
      bool SubbandAllocated = false;
      FlowToSchedule* scheduledFlow;
      int l_iScheduledFlowIndex = 0;

      for (int k = 0; k < flows->size (); k++)
        {
          if (metrics[s][k] > targetMetric && !l_bFlowScheduled[k])
            {
              targetMetric = metrics[s][k];
              SubbandAllocated = true;
              scheduledFlow = flows->at (k);
              l_iScheduledFlowIndex = k;
            }
        }

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
          //assert(scheduledFlow->m_bearer->GetApplication()->GetApplicationID() == l_iScheduledFlowIndex);
          int alloc_num = scheduledFlow->GetListOfAllocatedRBs()->size();
          int transportBlockSize = amc->GetTBSizeFromMCS (mcs, alloc_num);
          if (transportBlockSize >= scheduledFlow->GetDataToTransmit() * 8 )
          {
              //std::cout << "flow index: " << l_iScheduledFlowIndex << " alloc_rbs:" << alloc_num << std::endl;
              l_bFlowScheduled[l_iScheduledFlowIndex] = true;
              l_iScheduledFlows++;
          }
        }
    }

  delete [] l_bFlowScheduled;
  delete [] l_bFlowScheduledSINR;

  //Finalize the allocation
  PdcchMapIdealControlMessage *pdcchMsg = new PdcchMapIdealControlMessage ();

  for (FlowsToSchedule::iterator it = flows->begin (); it != flows->end (); it++)
    {
      FlowToSchedule *flow = (*it);

#ifdef SCHEDULER_DEBUG
	    std::cout << "Flow: " << flow->GetBearer()->GetApplication()->GetApplicationID();
	    for (int rb = 0; rb < flow->GetListOfAllocatedRBs()->size(); rb++) {
        int rbid = flow->GetListOfAllocatedRBs()->at(rb);
        if (rbid % rbg_size == 0)
          std::cout << " " << rbid / rbg_size;
	    }
	    std::cout << std::endl;
#endif

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
          //int transportBlockSize = amc->GetTBSizeFromMCS (mcs);
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
DownlinkNVSScheduler::ComputeSchedulingMetric(RadioBearer *bearer, double spectralEfficiency, int subChannel, double wbEff)
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
    default:
      metric = spectralEfficiency;
  }
  return metric;
}

void
DownlinkNVSScheduler::UpdateAverageTransmissionRate (int slice_serve)
{
  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
    {
	    RadioBearer *bearer = (*it);
      int app_id = bearer->GetApplication()->GetApplicationID();
      if (type2_app_[app_id] != slice_serve)
        continue;
      bearer->UpdateAverageTransmissionRate ();
    }
}
