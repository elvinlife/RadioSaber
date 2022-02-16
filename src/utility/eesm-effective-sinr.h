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


#ifndef EESM_EFFECTIVE_SINR_H_
#define EESM_EFFECTIVE_SINR_H_

#include <math.h>
#include <vector>

static double beta_value [20] = {
1.49,
1.53,
1.57,
1.61,
1.69,
1.69,
1.65,
3.36,
4.56,
6.42,
7.33,
7.68,
9.21,
10.81,
13.76,
17.52,
20.57,
22.75,
25.16,
28.38
};

static double
GetEesmEffectiveSinr (std::vector <double> &sinr)
{
  double eff_sinr;
  double sum_I_sinr = 0;
  double beta = 1;
  std::vector<double>::iterator it;
  for (it = sinr.begin (); it != sinr.end (); it++)
  {
	  //since sinr[] is expressed in dB we should convert it in natural unit!
	  double s = pow (10, ((*it)/10));
	  sum_I_sinr += exp (-s / beta);
  }

  eff_sinr = - beta * log (sum_I_sinr / sinr.size ());
  eff_sinr = 10 * log10 (eff_sinr); //convert in dB
  return eff_sinr;

  //double eff_sinr;
  //double sum_I_sinr = 0;
  //std::vector<double>::iterator it;
  //for (it = sinr.begin (); it != sinr.end (); it++)
  //  {
	//  //since sinr[] is expressed in dB we should convert it in natural unit!
	//  double s = pow (10, ((*it)/10));
	//  sum_I_sinr += s;
  //  }

  //eff_sinr = 10 * log10 (sum_I_sinr / sinr.size()); //convert in dB

  return eff_sinr;
}

static int
get_subband_size(int nof_prb)
{
  if (nof_prb <= 26) {
    return 4;
  }
  else if (nof_prb <= 63) {
    return 6;
  }
  else if (nof_prb <= 110) {
    return 8;
  }
  // 5G 100MHz
  // If a RB is 400khz, there're 250 PRB. Thus, the subband size is 16 PRB(32 200khz PRB)
  // If a RB is 800khz, there're 125 PRB. Thus, the subband size is 8 PRB(32 200khz PRB)
  else if (nof_prb <= 500) {
    return 32;
  }
  else {
    throw std::runtime_error("Num of physical RBs above 110");
  }
}

static int
get_rbg_size(int nof_prb)
{
  if (nof_prb <= 10)
    return 1;
  else if (nof_prb <= 26)
    return 2;
  else if (nof_prb <= 63)
    return 3;
  else if (nof_prb <= 110)
    return 4;
  // According to 3GPP TS 38.214, the RBG size can be 2 or 4 or wideband(which means cannot allocate non-contiguous RBs)
  // https://www.etsi.org/deliver/etsi_ts/136200_136299/136213/15.10.00_60/ts_136213v151000p.pdf
  else if (nof_prb <= 500)
    return 8;
  else {
    throw std::runtime_error("Num of physical RBs above 110");
  }
}

#endif /* EESM_EFFECTIVE_SINR_H_ */
