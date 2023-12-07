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

#ifndef RADIOBEARER_H_
#define RADIOBEARER_H_

#include <unordered_map>

#include "radio-bearer-instance.h"

class NetworkNode;
class ClassifierParameters;
class Application;
class MacQueue;
class QoSParameters;
class Packet;
class RlcEntity;

class RadioBearer : public RadioBearerInstance {
 public:
  RadioBearer();
  virtual ~RadioBearer();

  MacQueue* GetMacQueue(void);

  void SetApplication(Application* a);
  Application* GetApplication(void);
  int GetUserID(void);
  // higher value means higher priority
  int GetPriority(void);

  void UpdateTransmittedBytes(int bytes);
  int GetTransmittedBytes(void) const;
  void ResetTransmittedBytes(void);
  void UpdateAverageTransmissionRate();
  double GetAverageTransmissionRate(void) const;
  void SetLastUpdate(void);
  double GetLastUpdate(void) const;
  void UpdateCumulateRBs(int rbs);
  unsigned long GetCumulateRBs(void) const;
  unsigned long GetCumulateBytes(void) const;

  void Enqueue(Packet* packet);
  bool HasPackets(void);

  Packet* CreatePacket(int bytes);

  void CheckForDropPackets();

  int GetQueueSize(void);
  double GetHeadOfLinePacketDelay(void);
  int GetHeadOfLinePacketSize(void);
  int GetByte(int byte);  // for FLS scheduler

  std::unordered_map<int, double>& GetFlowEnqueueInfo();

 private:
  Application* m_application;

  MacQueue* m_macQueue;

  // Scheduler Information
  double m_averageTransmissionRate;  // unit: bps
  int m_transmittedBytes;
  unsigned long m_cumulativeRBs;
  unsigned long m_cumulativeBytes;
  double m_lastUpdate;

  std::unordered_map<int, double> m_flow_enqueueInfo;
};

#endif /* RADIOBEARER_H_ */
