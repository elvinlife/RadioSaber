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



#include "enb-mac-entity.h"
#include "../packet/Packet.h"
#include "../packet/packet-burst.h"
#include "AMCModule.h"
#include "../../core/idealMessages/ideal-control-messages.h"
#include "../../device/NetworkNode.h"
#include "packet-scheduler/packet-scheduler.h"
#include "../../device/UserEquipment.h"
#include "../../device/ENodeB.h"
#include "../../load-parameters.h"
#include "../../core/eventScheduler/simulator.h"
#include <fstream>
#include <cassert>
#include <sstream>
#define CQI_INTERVAL 40
#define MAX_UE_TRACE 158
#define MAX_TTI_TRACE 475

EnbMacEntity::EnbMacEntity ()
{
  SetAmcModule (new AMCModule ());
  SetDevice (NULL);
  m_downlinkScheduler = NULL;
  m_uplinkScheduler = NULL;
  #ifdef USE_REAL_TRACE
  std::string fname = path + "cqi-traces-noise0/mapping.config";
  std::ifstream ifs(fname, std::ifstream::in);
  int uid, tid;
  while (ifs >> uid >> tid) {
    m_userMapping.push_back(tid);
  }
  #endif
}


EnbMacEntity::~EnbMacEntity ()
{
  delete m_downlinkScheduler;
  delete m_uplinkScheduler;
  Destroy ();
}


void
EnbMacEntity::SetUplinkPacketScheduler (PacketScheduler* s)
{
  m_uplinkScheduler = s;
}


void
EnbMacEntity::SetDownlinkPacketScheduler (PacketScheduler* s)
{
  m_downlinkScheduler = s;
}


PacketScheduler*
EnbMacEntity::GetUplinkPacketScheduler (void)
{
  return m_uplinkScheduler;
}


PacketScheduler*
EnbMacEntity::GetDownlinkPacketScheduler (void)
{
  return m_downlinkScheduler;
}


void
EnbMacEntity::ReceiveCqiIdealControlMessage  (CqiIdealControlMessage* msg)
{
  #ifndef USE_REAL_TRACE
#ifdef TEST_CQI_FEEDBACKS
  std::cout << "ReceiveIdealControlMessage (MAC) from  " << msg->GetSourceDevice ()->GetIDNetworkNode ()
		  << " to " << msg->GetDestinationDevice ()->GetIDNetworkNode () << std::endl;
#endif

  CqiIdealControlMessage::CqiFeedbacks *cqi = msg->GetMessage ();

  UserEquipment* ue = (UserEquipment*) msg->GetSourceDevice ();
  ENodeB* enb = (ENodeB*) GetDevice ();
  ENodeB::UserEquipmentRecord* record = enb->GetUserEquipmentRecord (ue->GetIDNetworkNode ());

  if (record != NULL)
    {
      std::vector<int> cqiFeedback;
      for (CqiIdealControlMessage::CqiFeedbacks::iterator it = cqi->begin (); it != cqi->end (); it++)
        {
	      cqiFeedback.push_back ((*it).m_cqi);
        }

#ifdef TEST_CQI_FEEDBACKS
      std::cout << "\t CQI: ";
      for (int i = 0; i < cqiFeedback.size (); i++)
        {
	      std::cout << cqiFeedback.at (i) << " ";
        }
      std::cout << std::endl;
#endif

#ifdef AMC_MAPPING
      std::cout << "\t CQI: ";
      for (int i = 0; i < cqiFeedback.size (); i++)
        {
	      std::cout << cqiFeedback.at (i) << " ";
        }
      std::cout << std::endl;

      std::cout << "\t MCS: ";
      for (int i = 0; i < cqiFeedback.size (); i++)
        {
	      std::cout << GetAmcModule ()->GetMCSFromCQI (cqiFeedback.at (i)) << " ";
        }
      std::cout << std::endl;

      std::cout << "\t TB: ";
      for (int i = 0; i < cqiFeedback.size (); i++)
        {
	      std::cout << GetAmcModule ()->GetTBSizeFromMCS(
	    		  GetAmcModule ()->GetMCSFromCQI (cqiFeedback.at (i))) << " ";
        }
      std::cout << std::endl;
#endif

      record->SetCQI (cqiFeedback);

    }
  else
    {
      std::cout << "ERROR: received cqi from unknow ue!"<< std::endl;
    }
    #endif

    #ifdef USE_REAL_TRACE  
    CqiIdealControlMessage::CqiFeedbacks *cqi = msg->GetMessage ();
    int nb_rbs = cqi->size();
    // assert(nb_rbs == 500);
    UserEquipment* ue = (UserEquipment*) msg->GetSourceDevice ();
    int user_id = ue->GetIDNetworkNode();
    ENodeB* enb = (ENodeB*) GetDevice ();
    ENodeB::UserEquipmentRecord* record = enb->GetUserEquipmentRecord(user_id);

    if (m_userCQITrace.find(user_id) == m_userCQITrace.end()) {
      int trace_id = m_userMapping[user_id % m_userMapping.size()];
      std::cerr << "user " << user_id << " uses trace " << trace_id << std::endl;
      std::string fname = path + "cqi-traces-noise0/ue" + std::to_string(trace_id) + ".log";
      std::ifstream ifs(fname, std::ifstream::in);
      int cqi = 0;
      for (int i = 0; i < MAX_TTI_TRACE; ++i) {
        std::string line;
        std::getline(ifs, line);
        std::istringstream iss(line);
        std::vector<int> cqi_one_tti;
        for (int j = 0; j < nb_rbs; ++j) {
          iss >> cqi;
          cqi_one_tti.push_back(cqi);
        }
        m_userCQITrace[user_id].push_back(cqi_one_tti);
      }
      ifs.close();
    }

    auto& cqi_traces = m_userCQITrace.at(user_id);
    int time_stamp = Simulator::Init ()->Now ()*1000 / CQI_INTERVAL;
    record->SetCQI( cqi_traces[time_stamp % cqi_traces.size() ] );

    #endif
}


void
EnbMacEntity::SendPdcchMapIdealControlMessage  (PdcchMapIdealControlMessage* msg)
{
}


void
EnbMacEntity::ReceiveSchedulingRequestIdealControlMessage (SchedulingRequestIdealControlMessage* msg)
{
  UserEquipment* ue = (UserEquipment*) msg->GetSourceDevice ();
  ENodeB* enb = (ENodeB*) GetDevice ();
  ENodeB::UserEquipmentRecord* record = enb->GetUserEquipmentRecord (ue->GetIDNetworkNode ());

  int bufferStatusReport = msg->GetBufferStatusReport ();

  if (record != NULL)
	{
	  record->SetSchedulingRequest (bufferStatusReport);
	}
  else
    {
      std::cout << "ERROR: received cqi from unknow ue!"<< std::endl;
    }
}
