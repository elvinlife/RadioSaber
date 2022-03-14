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

#include "fullband-cqi-manager.h"
#include "../NetworkNode.h"
#include "../../core/idealMessages/ideal-control-messages.h"
#include "../../protocolStack/mac/AMCModule.h"
#include "../UserEquipment.h"
#include "../ENodeB.h"
#include "../HeNodeB.h"
#include "../../phy/lte-phy.h"
#include "../../core/spectrum/bandwidth-manager.h"
#include "../../load-parameters.h"
#include <exception>
#include <cassert>
#include "../../utility/eesm-effective-sinr.h"


FullbandCqiManager::FullbandCqiManager()
{}

FullbandCqiManager::~FullbandCqiManager()
{}



std::vector<double>
FullbandCqiManager::GenSubbandSINR( std::vector<double> sinr )
{
  int subband_size = get_subband_size(sinr.size());
  std::vector<double> subbands_sinr(sinr.size(), 0);
  int l = 0, r = 0;
  while (l < sinr.size()) {
    r = l + subband_size;
    if (r >= sinr.size()) r = sinr.size();
    std::vector<double> subset(sinr.cbegin() + l, sinr.cbegin() + r);
    //assert(subset.size() == subband_size);
    double effective_sinr = GetEesmEffectiveSinr(subset);
    for (int i = l; i < r; i++) {
      subbands_sinr[i] = effective_sinr;
    }
    l = r;
  }
  return subbands_sinr;
}

void
FullbandCqiManager::CreateCqiFeedbacks (std::vector<double> sinr)
{
#ifdef TEST_CQI_FEEDBACKS
  std::cout << "FullbandCqiManager -> CreateCqiFeedbacks " << std::endl;
#endif

#ifdef AMC_MAPPING
  std::cout << "\t sinr: ";
  for (int i = 0; i < sinr.size (); i++)
    {
	  std::cout << sinr.at (i) << " ";
    }
  std::cout << std::endl;
#endif

  // subband cqi feedbacks generated here
  UserEquipment* thisNode = (UserEquipment*) GetDevice ();
  NetworkNode* targetNode = thisNode->GetTargetNode ();

  AMCModule *amc = GetDevice ()->GetProtocolStack ()->GetMacEntity ()->GetAmcModule ();

  sinr = GenSubbandSINR(sinr);

  std::vector<int> cqi = amc->CreateCqiFeedbacks (sinr);

  //int subband_size = get_subband_size(sinr.size());
  //for (int i = 0; i < cqi.size();) {
  //  std::cerr << cqi.at(i) << " ";
  //  i += subband_size;
  //}
  //std:cerr << std::endl;

  CqiIdealControlMessage *msg = new CqiIdealControlMessage ();
  msg->SetSourceDevice (thisNode);
  msg->SetDestinationDevice (targetNode);

  int nbSubChannels = cqi.size ();
  std::vector<double> dlSubChannels = thisNode->GetPhy ()->GetBandwidthManager ()->GetDlSubChannels ();

  for (int i = 0; i < nbSubChannels; i++)
    {
      msg->AddNewRecord (dlSubChannels.at (i), cqi.at (i));
    }

  SetLastSent ();

  thisNode->GetPhy ()->SendIdealControlMessage (msg);
}
