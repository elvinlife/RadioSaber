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

#ifndef PACKETSCHEDULER_H_
#define PACKETSCHEDULER_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "../../../core/idealMessages/ideal-control-messages.h"

const int MAX_BEARERS = 2;

class MacEntity;
class PacketBurst;
class Packet;
class RadioBearer;

struct SchedulerAlgoParam {
  int alpha;
  int beta;
  int epsilon;
  int psi;
  SchedulerAlgoParam(int _alpha, int _beta, int _epsilon, int _psi) {
    alpha = _alpha;
    beta = _beta;
    epsilon = _epsilon;
    psi = _psi;
  }
};

class PacketScheduler {
 public:
  struct FlowToSchedule {
    FlowToSchedule(RadioBearer* bearer, int dataToTransmit);
    virtual ~FlowToSchedule();
    RadioBearer* m_bearer;
    int m_allocatedBits;    // bits
    int m_transmittedData;  // bytes
    int m_dataToTransmit;   // bytes

    std::vector<double> m_spectralEfficiency;
    std::vector<int> m_listOfAllocatedRBs;
    std::vector<int> m_listOfSelectedMCS;
    std::vector<int> m_cqiFeedbacks;
    int m_wideBandCQI;

    RadioBearer* GetBearer(void);

    void UpdateAllocatedBits(int allocatedBits);
    int GetAllocatedBits(void) const;
    int GetTransmittedData(void) const;
    void SetDataToTransmit(int dataToTransmit);
    int GetDataToTransmit(void) const;

    void SetSpectralEfficiency(std::vector<double>& s);
    std::vector<double> GetSpectralEfficiency(void);

    void SetWidebandCQI(int);
    int GetWidebandCQI(void);

    std::vector<int>* GetListOfAllocatedRBs();
    std::vector<int>* GetListOfSelectedMCS();

    void SetCqiFeedbacks(std::vector<int>& cqiFeedbacks);
    std::vector<int> GetCqiFeedbacks(void);
  };

  struct UserToSchedule {
   private:
    int m_userID;
    NetworkNode* m_userNode;
    int m_allocatedBits;

    // spectralEfficiency is transmission rate per HZ(bps/HZ)
    std::vector<double> m_spectralEfficiency;
    std::vector<int> m_listOfAllocatedRBs;
    int m_wideBandCQI;
    std::vector<int> m_cqiFeedbacks;

   public:
    RadioBearer* m_bearers[MAX_BEARERS];
    int m_dataToTransmit[MAX_BEARERS];
    int m_requiredRBs;
    UserToSchedule(int, NetworkNode*);
    virtual ~UserToSchedule();

    int GetUserID(void) { return m_userID; }
    NetworkNode* GetUserNode(void) { return m_userNode; }
    void UpdateAllocatedBits(int allocatedBits);
    int GetAllocatedBits(void);
    std::vector<RadioBearer*> GetBearers(void);

    void SetSpectralEfficiency(std::vector<double>& s);
    std::vector<double>& GetSpectralEfficiency(void);
    void SetWidebandCQI(int);
    int GetWidebandCQI(void);
    int GetRequiredRBs(void);
    void SetCqiFeedbacks(std::vector<int>& cqiFeedbacks);
    std::vector<int>& GetCqiFeedbacks(void);

    std::vector<int>* GetListOfAllocatedRBs();
    double GetAverageTransmissionRate();
  };

  PacketScheduler();
  virtual ~PacketScheduler();

  void Destroy(void);

  void SetMacEntity(MacEntity* mac);
  MacEntity* GetMacEntity(void);

  void Schedule(void);
  virtual void DoSchedule(void);

  void StopSchedule();
  virtual void DoStopSchedule();

  void InsertFlowToSchedule(RadioBearer* bearer, int dataToTransmit,
                            std::vector<double> specEff,
                            std::vector<int> cqiFeedbacks);

  void UpdateAllocatedBits(FlowToSchedule* scheduledFlow, int allocatedBits,
                           int allocatedRB, int selectedMCS);

  void CheckForDLDropPackets();

  void UpdateTimeStamp();

  unsigned long GetTimeStamp();

  typedef std::vector<FlowToSchedule*> FlowsToSchedule;
  void CreateFlowsToSchedule(void);
  void DeleteFlowsToSchedule(void);
  void ClearFlowsToSchedule();
  FlowsToSchedule* GetFlowsToSchedule(void) const;

  // user-aware scheduling
  // typedef std::unordered_map<int, UserToSchedule*> UsersToSchedule;
  typedef std::vector<UserToSchedule*> UsersToSchedule;
  void InsertFlowToUser(RadioBearer* bearer, int dataToTransmit,
                        std::vector<double> specEff,
                        std::vector<int> TestCqiFeedbacks);
  void CreateUsersToSchedule(void);
  void DeleteUsersToSchedule(void);
  void ClearUsersToSchedule();
  UsersToSchedule* GetUsersToSchedule(void) const;

 private:
  MacEntity* m_mac;
  FlowsToSchedule* m_flowsToSchedule;
  UsersToSchedule* m_usersToSchedule;
  unsigned long m_ts;
};

#endif /* PACKETSCHEDULER_H_ */
