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


#include "downlink-transport-scheduler.h"
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
#include <fstream>
#include <cassert>
#include <algorithm>
using std::vector;

using coord_t = std::pair<int, int>;
using coord_cqi_t = std::pair<coord_t, int>;

DownlinkTransportScheduler::DownlinkTransportScheduler(std::string config_fname, int algo)
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
  inter_sched_ = algo;
  SetMacEntity (0);
  CreateFlowsToSchedule ();
}

DownlinkTransportScheduler::~DownlinkTransportScheduler()
{
  Destroy ();
}

void DownlinkTransportScheduler::SelectFlowsToSchedule ()
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
		  InsertFlowToSchedule(bearer, dataToTransmit, spectralEfficiency, cqiFeedbacks);
		}
	  else
	    {}
	}
  //std::cerr << std::endl;
}

void
DownlinkTransportScheduler::DoSchedule (void)
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
DownlinkTransportScheduler::DoStopSchedule (void)
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

static vector<int> MaximizeCell(int** flow_cqi, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  // it's ok that the quota is negative
  vector<coord_cqi_t> sorted_cqi;
  vector<int>         slice_rbgs(nb_slices, 0);
  vector<int>         rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i)
    for (int j = 0; j < nb_slices; ++j) {
      sorted_cqi.emplace_back(coord_t(i, j), flow_cqi[i][j]);
    }
  std::sort(sorted_cqi.begin(), sorted_cqi.end(), [](coord_cqi_t a, coord_cqi_t b){return a.second > b.second; });
  for (auto it = sorted_cqi.begin(); it != sorted_cqi.end(); ++it) {
    int rbg_id = it->first.first;
    int slice_id = it->first.second;
    if (slice_rbgs[slice_id] < slice_quota_rbgs[slice_id] && rbg_to_slice[rbg_id] == -1) {
      rbg_to_slice[rbg_id] = slice_id;
      slice_rbgs[slice_id] += 1;
    } 
  }
  return rbg_to_slice;
}

static vector<int> VogelApproximate(int** flow_cqi, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  vector<int>   slice_rbgs(nb_slices, 0);
  vector<int>   rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i) {
    int max_diff = -1;
    coord_t coord_1st;
    coord_t coord_2nd;
    // horizonal search
    for (int j = 0; j < nb_rbgs; ++j) {
      // the rbg has been allocated
      if (rbg_to_slice[j] != -1)
        continue;
      int cqi_1st = -1, cqi_2nd = -1;
      int slice_1st, slice_2nd;

      for (int k = 0; k < nb_slices; ++k) {
        // the slice has reached quota
        if (slice_rbgs[k] >= slice_quota_rbgs[k])
          continue;
        if (cqi_1st == -1 || flow_cqi[j][k] > cqi_1st) {
          slice_1st = k;
          cqi_1st = flow_cqi[j][k];
          continue;
        }
        if (cqi_2nd == -1 || flow_cqi[j][k] > cqi_2nd) {
          slice_2nd = k;
          cqi_2nd = flow_cqi[j][k];
          continue;
        }
      }
      if (cqi_1st - cqi_2nd > max_diff) {
        max_diff = cqi_1st - cqi_2nd;
        coord_1st = coord_t(j, slice_1st);
        coord_2nd = coord_t(j, slice_2nd);
      }
    }
    // vertical search
    for (int k = 0; k < nb_slices; ++k) {
      // the slice has reached quota
      if (slice_rbgs[k] >= slice_quota_rbgs[k])
        continue;
      int cqi_1st = -1, cqi_2nd = -1;
      int rbg_1st, rbg_2nd;
      for (int j = 0; j < nb_rbgs; ++j) {
        if (rbg_to_slice[j] != -1)
          continue;
        if (cqi_1st == -1 || flow_cqi[j][k] > cqi_1st) {
          rbg_1st = j;
          cqi_1st = flow_cqi[j][k];
          continue;
        }
        if (cqi_2nd == -1 || flow_cqi[j][k] > cqi_2nd) {
          rbg_2nd = j;
          cqi_2nd = flow_cqi[j][k];
          continue;
        }
      }
      if (cqi_1st - cqi_2nd > max_diff) {
        max_diff = cqi_1st - cqi_2nd;
        coord_1st = coord_t(rbg_1st, k);
        coord_2nd = coord_t(rbg_2nd, k);
      }
    }
    rbg_to_slice[coord_1st.first] = coord_1st.second;
    slice_rbgs[coord_1st.second] += 1;
  }
  return rbg_to_slice;
}

void
DownlinkTransportScheduler::RBsAllocation ()
{
#ifdef SCHEDULER_DEBUG
	std::cout << " ---- DownlinkTransportScheduler::RBsAllocation";
#endif
  FlowsToSchedule* flows = GetFlowsToSchedule ();
  int nb_rbs = GetMacEntity ()->GetDevice ()->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ().size ();
  int rbg_size = get_rbg_size(nb_rbs);
  // we only evaluate 20Mhz with 100rbs here
  assert(rbg_size == 4 && nb_rbs == 100);
  int nb_rbgs = (nb_rbs + rbg_size - 1) / rbg_size;

  std::vector<double> slice_target_rbs(num_type2_slices_);
  std::vector<int> slice_final_rbgs(num_type2_slices_, 0);
  std::vector<int> slice_quota_rbgs(num_type2_slices_, 0);
  int extra_rbgs = nb_rbgs;
  std::cout << "\t slice target RBs:";
  for (int i = 0; i < num_type2_slices_; ++i) {
    double quota = nb_rbs * type2_weights_[i];
    slice_target_rbs[i] = quota + type2_rbs_offset_[i];
    slice_quota_rbgs[i] = (int)(slice_target_rbs[i] / rbg_size);
    extra_rbgs -= slice_quota_rbgs[i];
    std::cout << "\t" << slice_target_rbs[i];
  }
  std::cout << std::endl;
  while (extra_rbgs > 0) {
    slice_quota_rbgs[rand() % num_type2_slices_] += 1;
    extra_rbgs -= 1;
  }

  // create a matrix of flow metrics (RBG, flow index)
  double metrics[nb_rbgs][flows->size ()];

  for (int i = 0; i < nb_rbgs; i++) {
	  for (int j = 0; j < flows->size (); j++) {
		  metrics[i][j] = ComputeSchedulingMetric (
        flows->at (j)->GetBearer (),
        flows->at (j)->GetSpectralEfficiency ().at (i * rbg_size),
        i, flows->at(j)->GetAllEfficiency());
	  }
  }

#ifdef SCHEDULER_DEBUG
  //std::cout << ", available RBGs " << nb_rbgs << ", flows " << flows->size () << std::endl;
  for (int ii = 0; ii < flows->size (); ii++)
  {
	std::cout << "\t metrics for flow "
		  << flows->at (ii)->GetBearer ()->GetApplication ()->GetApplicationID () << ":";
	for (int jj = 0; jj < nb_rbgs; jj++) {
      fprintf(stdout, " (%d, %.3f, %d, %.3f)",
          jj, metrics[jj][ii], 
          flows->at(ii)->GetCqiFeedbacks().at(jj * rbg_size),
          flows->at(ii)->GetSpectralEfficiency().at(jj * rbg_size));
	  }
	std::cout << std::endl;
  }
#endif

  AMCModule *amc = GetMacEntity ()->GetAmcModule ();

  // Assumption: every slice has at least one flow
  int **flow_id = new int*[nb_rbgs];
  int **flow_cqi = new int*[nb_rbgs];

  // Assign the flow_id and cqi to every slice in every rbg
  for (int s = 0; s < nb_rbgs; s++)
  {
    flow_id[s] = new int[num_type2_slices_];
    flow_cqi[s] = new int[num_type2_slices_];
    std::vector<double> max_ranks(num_type2_slices_, -1);     // the highest flow metric in every slice
    for (int k = 0; k < flows->size(); k++)
    {
      int app_id = flows->at(k)->GetBearer()->GetApplication()->GetApplicationID();
      int slice_id = type2_app_[app_id];
      if (metrics[s][k] > max_ranks[slice_id]) {
        max_ranks[slice_id] = metrics[s][k];
        flow_id[s][slice_id] = k;
        flow_cqi[s][slice_id] = flows->at(k)->GetSpectralEfficiency().at(s * rbg_size);
      }
    }
  }

  vector<int> rbg_to_slice;
  if (inter_sched_ == 1) {
    rbg_to_slice = VogelApproximate(flow_cqi, slice_quota_rbgs, nb_rbgs, num_type2_slices_);
  }
  else {
    rbg_to_slice = MaximizeCell(flow_cqi, slice_quota_rbgs, nb_rbgs, num_type2_slices_);
  }

  for (int i = 0; i < rbg_to_slice.size(); ++i) {
    int fid = flow_id[i][rbg_to_slice[i]];
    int sid = type2_app_[flows->at(fid)->GetBearer()->GetApplication()->GetApplicationID()];
    slice_final_rbgs[sid] += 1;
    int l = i * rbg_size, r = (i+1) * rbg_size;
    for (int j = l; j < r; ++j) {
      flows->at(fid)->GetListOfAllocatedRBs()->push_back(j);
    }
  }

  for (int i = 0; i < num_type2_slices_; ++i) {
    type2_rbs_offset_[i] = slice_target_rbs[i] - slice_final_rbgs[i] * rbg_size;
  }


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
DownlinkTransportScheduler::ComputeSchedulingMetric(RadioBearer *bearer, double spectralEfficiency, int subChannel, double wbEff)
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
DownlinkTransportScheduler::UpdateAverageTransmissionRate (void)
{
  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
    {
      RadioBearer *bearer = (*it);
      bearer->UpdateAverageTransmissionRate ();
    }
}