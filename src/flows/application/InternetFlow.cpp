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

#include "InternetFlow.h"
#include <cstdlib>
#include "../../componentManagers/NetworkManager.h"
#include "../radio-bearer.h"
#include <cmath>

// 1460.000000,0.500000
// 2920.000000,0.600000
// 4380.000000,0.700000
// 7300.000000,0.750000
// 10220.000000,0.800000
// 58400.000000, 0.812500
// 105120.000000, 0.825000
// 200020.000000, 0.850000
// 389820.000000,0.900000
// 1733020.000000, 0.95000
// 3076220.000000,1.000000

#define MAXMTUSIZE 1490
int get_avg_flowsize() {
  double avg_flowsize = 0;
  for (int i = 0; i < InternetFlow::m_typeflow; i++) {
    if (i == 0)
      avg_flowsize += InternetFlow::m_flowcdf[i] * InternetFlow::m_flowsize[i];
    else
      avg_flowsize += (InternetFlow::m_flowcdf[i] - InternetFlow::m_flowcdf[i-1]) * InternetFlow::m_flowsize[i];
  }
  return (int)avg_flowsize;
}

int InternetFlow::m_typeflow = 11;
int InternetFlow::m_flowsize[] = {1460, 2920, 4380, 7300, 10220, 58400, 105120, 200020, 389820, 1733020, 3076220};
double InternetFlow::m_flowcdf[] = {0.5, 0.6, 0.7, 0.75, 0.8, 0.8125, 0.825, 0.85, 0.9, 0.95, 1};
int InternetFlow::m_avg_flowsize = get_avg_flowsize();

InternetFlow::InternetFlow()
{
  SetApplicationType (Application::APPLICATION_TYPE_IPFLOW);
  m_flowCounter = 0;
}

InternetFlow::~InternetFlow()
{
  Destroy ();
}

void
InternetFlow::DoStart (void)
{
  Simulator::Init()->Schedule(0.0, &InternetFlow::Send, this);
}

void
InternetFlow::DoStop (void)
{
}

void
InternetFlow::ScheduleTransmit (double time)
{
  if ( (Simulator::Init()->Now () + time) < GetStopTime () )
    {
      Simulator::Init()->Schedule(time, &InternetFlow::Send, this);
    }
}

void
InternetFlow::Send (void)
{
  //CREATE A NEW PACKET (ADDING UDP, IP and PDCP HEADERS)
  int flow_size = GetSize();
  uint16_t last_pkt_size = flow_size % MAXMTUSIZE;
  int n_pkts = std::ceil( flow_size / (double)MAXMTUSIZE );
  vector<Packet*> packets;
  for (int i = 0; i < n_pkts; i++) {
    Packet *packet = new Packet ();
    int uid = Simulator::Init()->GetUID ();
    packet->SetID(uid);
    packet->SetTimeStamp (Simulator::Init()->Now ());
    packet->SetSize (MAXMTUSIZE);
    PacketTAGs *tags = new PacketTAGs();
    tags->SetApplicationType(PacketTAGs::APPLICATION_TYPE_IPFLOW);
    tags->SetFrameNumber(m_flowCounter);
    packet->SetPacketTags(tags);

    UDPHeader* udp = new UDPHeader(
      GetClassifierParameters()->GetSourcePort(),
      GetClassifierParameters()->GetDestinationPort()
      );
    packet->AddUDPHeader(udp);
    IPHeader* ip = new IPHeader(
      GetClassifierParameters()->GetSourceID(),
      GetClassifierParameters()->GetDestinationID()
    );
    packet->AddIPHeader(ip);
    PDCPHeader* pdcp = new PDCPHeader();
    packet->AddPDCPHeader(pdcp);

    packets.push_back(packet);
  }

  // set the correct pkt size of the last pkt
  if (last_pkt_size != 0) {
    packets.back()->SetSize(last_pkt_size);
  }
  packets.front()->GetPacketTags()->SetStartByte(1);
  packets.front()->GetPacketTags()->SetApplicationSize(flow_size);
  packets.back()->GetPacketTags()->SetEndByte(1);
  packets.back()->GetPacketTags()->SetApplicationSize(flow_size);
  for (auto it = packets.begin(); it != packets.end(); ++it) {
    GetRadioBearer()->Enqueue(*it);
  }

  m_flowCounter += 1;
  ScheduleTransmit( GetInterval() );
}

void
InternetFlow::SetAvgRate(double rate)
{
  // Bytes / Mbps => *8 (us) /1000000 (s)
  m_interval = InternetFlow::m_avg_flowsize / rate * 8 / 1000000;
  m_lambda = 1 / m_interval;
  m_distribute = std::exponential_distribution<double>(m_lambda);
  // std::cerr << "IPFlow rate: " << rate
  //     << " mbps; lambda: " << m_lambda
  //     << " ; flow size: " << InternetFlow::m_avg_flowsize << std::endl;
}

double
InternetFlow::GetInterval(void)
{
  // return m_interval;
  double interval = m_distribute(m_generator);
  return std::ceil(interval * 1000) / 1000.0;
}

int
InternetFlow::GetSize(void) const
{
  double cdf = (double) rand() / RAND_MAX;
  for (int i = 0; i < InternetFlow::m_typeflow; i++) {
    if (InternetFlow::m_flowcdf[i] >= cdf) {
      return InternetFlow::m_flowsize[i];
    }
  }
  throw std::runtime_error("CDF should be in [0, 1]");
}