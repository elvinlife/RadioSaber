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

#ifndef INTERNETFLOW_H_
#define INTERNETFLOW_H_

#include <random>

#include "Application.h"

class InternetFlow : public Application {
 public:
  static int m_avg_flowsize;
  static int m_typeflow;
  static int m_flowsize[];
  static double m_flowcdf[];

  InternetFlow();
  virtual ~InternetFlow();

  virtual void DoStart(void);
  virtual void DoStop(void);

  void ScheduleTransmit(double time);
  void Send(void);
  // set average sending rate(in Mbps)
  void SetAvgRate(double rate);

 private:
  double GetInterval(void);
  int GetSize(void) const;
  double m_interval;
  int m_flowCounter;

  double m_lambda;
  std::default_random_engine m_generator;
  std::exponential_distribution<double> m_distribute;
};

#endif /* CBR_H_ */
