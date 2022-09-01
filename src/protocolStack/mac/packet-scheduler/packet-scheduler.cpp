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


#include "packet-scheduler.h"
#include "../mac-entity.h"
#include "../../packet/Packet.h"
#include "../../packet/packet-burst.h"
#include "../../../device/NetworkNode.h"
#include "../../../flows/radio-bearer.h"
#include "../../../protocolStack/rrc/rrc-entity.h"
#include "../../../protocolStack/mac/AMCModule.h"
#include "../../../flows/application/Application.h"
#include "../../../flows/MacQueue.h"
#include "../../../flows/QoS/QoSParameters.h"
#include "../../rlc/am-rlc-entity.h"
#include "../../../utility/eesm-effective-sinr.h"
#include <cassert>

PacketScheduler::PacketScheduler()
{
  m_mac = NULL;
  m_flowsToSchedule = NULL;
  m_ts = 0;
}

PacketScheduler::~PacketScheduler()
{
  ClearFlowsToSchedule ();
  delete m_flowsToSchedule;
  m_mac = NULL;
}

void
PacketScheduler::Destroy (void)
{
  ClearFlowsToSchedule ();
  delete m_flowsToSchedule;
  m_mac = NULL;
}

void
PacketScheduler::SetMacEntity (MacEntity* mac)
{
  m_mac = mac;
}

MacEntity*
PacketScheduler::GetMacEntity (void)
{
  return m_mac;
}

void
PacketScheduler::Schedule (void)
{
  DoSchedule ();
}

void
PacketScheduler::DoSchedule (void)
{}

void
PacketScheduler::StopSchedule ()
{
  DoStopSchedule ();
}

void
PacketScheduler::DoStopSchedule ()
{}

PacketScheduler::FlowToSchedule::FlowToSchedule(RadioBearer* bearer, int dataToTransmit)
{
  m_bearer = bearer;
  m_allocatedBits = 0;
  m_transmittedData = 0;
  m_dataToTransmit = dataToTransmit;
}

PacketScheduler::FlowToSchedule::~FlowToSchedule()
{
}

void
PacketScheduler::CreateFlowsToSchedule (void)
{
  m_flowsToSchedule = new FlowsToSchedule ();
}

void
PacketScheduler::DeleteFlowsToSchedule (void)
{
  ClearFlowsToSchedule ();
  delete  m_flowsToSchedule;
}

PacketScheduler::FlowsToSchedule*
PacketScheduler::GetFlowsToSchedule (void) const
{
  return m_flowsToSchedule;
}

void
PacketScheduler::ClearFlowsToSchedule ()
{
  FlowsToSchedule*  records = GetFlowsToSchedule ();
  FlowsToSchedule::iterator iter;
  if (records == NULL)
    throw std::runtime_error("FLowstoSchedule is Null");

  for (iter = records->begin(); iter != records->end (); iter++)
  {
    delete *iter;
  }

  GetFlowsToSchedule ()->clear ();
}

RadioBearer*
PacketScheduler::FlowToSchedule::GetBearer (void)
{
  return m_bearer;
}

void
PacketScheduler::FlowToSchedule::SetSpectralEfficiency (std::vector<double>& s)
{
  m_spectralEfficiency = s;
}

std::vector<double>
PacketScheduler::FlowToSchedule::GetSpectralEfficiency (void)
{
  return m_spectralEfficiency;
}

void
PacketScheduler::FlowToSchedule::SetWidebandCQI(int cqi)
{
  m_wideBandCQI = cqi;
}

int
PacketScheduler::FlowToSchedule::GetWidebandCQI()
{
  return m_wideBandCQI;
}

void
PacketScheduler::FlowToSchedule::UpdateAllocatedBits (int allocatedBits)
{
  m_allocatedBits += allocatedBits;
  int availableBytes = m_allocatedBits/8;

  int transmittedPackets = ceil
			  (availableBytes/1513.0);

  // I guess this is to calculate the payload  
  m_transmittedData = availableBytes - (transmittedPackets * 8);
  if (m_transmittedData < 0)
  {
	  m_transmittedData = 0;
  }
}

int
PacketScheduler::FlowToSchedule::GetAllocatedBits (void) const
{
  return m_allocatedBits;
}

int
PacketScheduler::FlowToSchedule::GetTransmittedData (void) const
{
  return m_transmittedData;
}

void
PacketScheduler::FlowToSchedule::SetDataToTransmit (int dataToTransmit)
{
  m_dataToTransmit = dataToTransmit;
}

int
PacketScheduler::FlowToSchedule::GetDataToTransmit (void) const
{
  return m_dataToTransmit;
}

std::vector<int>*
PacketScheduler::FlowToSchedule::GetListOfAllocatedRBs ()
{
  return &m_listOfAllocatedRBs;
}

std::vector<int>*
PacketScheduler::FlowToSchedule::GetListOfSelectedMCS ()
{
  return &m_listOfSelectedMCS;
}

void
PacketScheduler::FlowToSchedule::SetCqiFeedbacks (std::vector<int>& cqiFeedbacks)
{
  m_cqiFeedbacks = cqiFeedbacks;
}

std::vector<int>
PacketScheduler::FlowToSchedule::GetCqiFeedbacks (void)
{
  return m_cqiFeedbacks;
}

void
PacketScheduler::InsertFlowToSchedule (RadioBearer* bearer, int dataToTransmit, std::vector<double> specEff, std::vector<int> cqiFeedbacks)
{
  FlowToSchedule *flowToSchedule = new FlowToSchedule(bearer, dataToTransmit);
  flowToSchedule->SetSpectralEfficiency (specEff);
  flowToSchedule->SetCqiFeedbacks (cqiFeedbacks);

  std::vector<double> sinrs;
	int numberOfCqi = cqiFeedbacks.size ();
  AMCModule *amc = GetMacEntity()->GetAmcModule();
	for (int i = 0; i < numberOfCqi; i++)
	{
    sinrs.push_back(amc->GetSinrFromCQI(cqiFeedbacks.at(i)));
	}
  double wideSINR = GetEesmEffectiveSinr(sinrs);
  flowToSchedule->SetWidebandCQI(amc->GetCQIFromSinr(wideSINR));

#ifdef SCHEDULER_DEBUG
	std::cout << "\t  --> selected flow: "
			<< bearer->GetApplication ()->GetApplicationID ()
			<< " data:" << dataToTransmit
      << " sinr:" << wideSINR << std::endl;
#endif

  GetFlowsToSchedule ()->push_back(flowToSchedule);
}

void
PacketScheduler::UpdateAllocatedBits (FlowToSchedule* scheduledFlow, int allocatedBits,  int allocatedRB, int selectedMCS)
{
  scheduledFlow->UpdateAllocatedBits (allocatedBits);
  scheduledFlow->GetListOfAllocatedRBs()->push_back(allocatedRB);
  scheduledFlow->GetListOfSelectedMCS()->push_back(selectedMCS);
}

void
PacketScheduler::CheckForDLDropPackets ()
{
  RrcEntity *rrc = GetMacEntity ()->GetDevice ()->GetProtocolStack ()->GetRrcEntity ();
  RrcEntity::RadioBearersContainer* bearers = rrc->GetRadioBearerContainer ();

  for (std::vector<RadioBearer* >::iterator it = bearers->begin (); it != bearers->end (); it++)
    {
	  //delete packets from queue
	  (*it)->GetMacQueue ()->CheckForDropPackets (
			  (*it)->GetQoSParameters ()->GetMaxDelay (), (*it)->GetApplication ()->GetApplicationID ());

	  //delete fragment waiting in AM RLC entity
	  if ((*it)->GetRlcEntity()->GetRlcModel() == RlcEntity::AM_RLC_MODE)
	    {
		  AmRlcEntity* amRlc = (AmRlcEntity*) (*it)->GetRlcEntity();
		  amRlc->CheckForDropPackets (
				  (*it)->GetQoSParameters ()->GetMaxDelay (), (*it)->GetApplication ()->GetApplicationID ());
	    }
    }
}

void
PacketScheduler::UpdateTimeStamp()
{
  m_ts += 1;
}

unsigned long
PacketScheduler::GetTimeStamp()
{
  return m_ts;
}

// dataToTransmit is in unit of bytes
void
PacketScheduler::InsertFlowToUser (RadioBearer* bearer, int dataToTransmit, std::vector<double> specEff, std::vector<int> cqiFeedbacks)
{
  int userID = bearer->GetUserID();
  int bearer_priority = bearer->GetPriority();

  for (auto it = m_usersToSchedule->begin(); it != m_usersToSchedule->end(); ++it) {
    if ((*it)->GetUserID() == userID) {
      (*it)->m_bearers[bearer_priority] = bearer;
      (*it)->m_dataToTransmit[bearer_priority] = dataToTransmit;
      return;
    }
  }
  UserToSchedule* user = new UserToSchedule(userID, bearer->GetDestination());
  user->SetSpectralEfficiency(specEff);
  user->SetCqiFeedbacks(cqiFeedbacks);
  
  std::vector<double> sinrs;
  AMCModule *amc = GetMacEntity()->GetAmcModule();
	for (int i = 0; i < cqiFeedbacks.size(); i++)
	{
    sinrs.push_back(amc->GetSinrFromCQI(cqiFeedbacks.at(i)));
	}
  int wideCQI = amc->GetCQIFromSinr(GetEesmEffectiveSinr(sinrs));
  user->SetWidebandCQI(wideCQI);

  assert(user->m_bearers[bearer_priority] == NULL);
  user->m_bearers[bearer_priority] = bearer;
  user->m_dataToTransmit[bearer_priority] = dataToTransmit;
  m_usersToSchedule->push_back(user);
  user->m_requiredRBs += (dataToTransmit * 8 / amc->GetTBSizeFromMCS(amc->GetMCSFromCQI(wideCQI))); 
}

void
PacketScheduler::CreateUsersToSchedule (void)
{
  m_usersToSchedule = new UsersToSchedule();
}

void
PacketScheduler::DeleteUsersToSchedule (void)
{
  ClearUsersToSchedule();
  delete m_usersToSchedule;
}

void
PacketScheduler::ClearUsersToSchedule()
{
  for (auto it = m_usersToSchedule->begin(); it != m_usersToSchedule->end(); ++it) {
    delete *it;
  }
  m_usersToSchedule->clear();
}

PacketScheduler::UsersToSchedule*
PacketScheduler::GetUsersToSchedule (void) const
{
  return m_usersToSchedule;
}

PacketScheduler::UserToSchedule::UserToSchedule(int id, NetworkNode* node)
{
  m_userID = id;
  m_userNode = node;
  m_allocatedBits = 0;
  m_wideBandCQI = 0;
  m_requiredRBs = 0;
  for (int i = 0; i < MAX_BEARERS; ++i) {
    m_bearers[i] = NULL;
    m_dataToTransmit[i] = 0;
  }
}

PacketScheduler::UserToSchedule::~UserToSchedule()
{}

void
PacketScheduler::UserToSchedule::SetSpectralEfficiency (std::vector<double>& s)
{
  m_spectralEfficiency = s;
}

std::vector<double>&
PacketScheduler::UserToSchedule::GetSpectralEfficiency (void)
{
  return m_spectralEfficiency;
}

void
PacketScheduler::UserToSchedule::SetCqiFeedbacks (std::vector<int>& cqiFeedbacks)
{
  m_cqiFeedbacks = cqiFeedbacks;
}

std::vector<int>&
PacketScheduler::UserToSchedule::GetCqiFeedbacks (void)
{
  return m_cqiFeedbacks;
}

void
PacketScheduler::UserToSchedule::SetWidebandCQI(int cqi)
{
  m_wideBandCQI = cqi;
}

int
PacketScheduler::UserToSchedule::GetWidebandCQI()
{
  return m_wideBandCQI;
}

std::vector<int>*
PacketScheduler::UserToSchedule::GetListOfAllocatedRBs ()
{
  return &m_listOfAllocatedRBs;
}

double
PacketScheduler::UserToSchedule::GetAverageTransmissionRate()
{
  double sum_rate = 1;
  for (int i = 0; i < MAX_BEARERS; i++) {
    if (m_bearers[i]) {
      sum_rate += m_bearers[i]->GetAverageTransmissionRate();
    }
  }
  return sum_rate;
}

void
PacketScheduler::UserToSchedule::UpdateAllocatedBits (int allocatedBits)
{
  m_allocatedBits += allocatedBits;
}

int
PacketScheduler::UserToSchedule::GetAllocatedBits()
{
  return m_allocatedBits;
}