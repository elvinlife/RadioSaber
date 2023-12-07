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

#ifndef MW_RULE_DOWNLINK_PACKET_SCHEDULER_H_
#define MW_RULE_DOWNLINK_PACKET_SCHEDULER_H_

#include "downlink-packet-scheduler.h"

class MwRulePacketScheduler : public DownlinkPacketScheduler {
 public:
  MwRulePacketScheduler();
  virtual ~MwRulePacketScheduler();

  virtual double ComputeSchedulingMetric(RadioBearer *bearer,
                                         double spectralEfficiency,
                                         int subChannel);
};

#endif /* MW_RULE_DOWNLINK_PACKET_SCHEDULER_H_ */
