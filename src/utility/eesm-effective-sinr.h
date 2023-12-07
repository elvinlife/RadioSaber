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

static double beta_value[20] = {1.49,  1.53,  1.57,  1.61,  1.69,  1.69, 1.65,
                                3.36,  4.56,  6.42,  7.33,  7.68,  9.21, 10.81,
                                13.76, 17.52, 20.57, 22.75, 25.16, 28.38};

static double GetEesmEffectiveSinr(std::vector<double> &sinr) {
  double eff_sinr;
  double sum_I_sinr = 0;
  double beta = 1;
  for (auto it = sinr.begin(); it != sinr.end(); it++) {
    // since sinr[] is expressed in dB we should convert it in natural unit!
    double s = pow(10, ((*it) / 10));
    sum_I_sinr += exp(-s / beta);
  }

  eff_sinr = -beta * log(sum_I_sinr / sinr.size());
  eff_sinr = 10 * log10(eff_sinr);  // convert in dB
  return eff_sinr;

  // double eff_sinr;
  // double sum_I_sinr = 0;
  // double beta = 1;
  // std::vector<double>::iterator it;
  // for (it = sinr.begin (); it != sinr.end (); it++)
  // {
  //   //use sinr in db
  //   double s = *it;
  //   sum_I_sinr += exp (-s / beta);
  // }

  // eff_sinr = - beta * log (sum_I_sinr / sinr.size ());
  // return eff_sinr;
}

static int get_subband_size(int nof_prb) {
  if (nof_prb <= 26) {
    return 4;
  } else if (nof_prb <= 63) {
    return 6;
  } else if (nof_prb <= 110) {
    return 8;
  }
  // 5G 100MHz
  // https://www.etsi.org/deliver/etsi_ts/136200_136299/136213/15.10.00_60/ts_136213v151000p.pdf
  // Table 7.2.1-3 64-110 RBs => Subband Size(8 RBs) If the subcarrier is 60khz,
  // a RB is 800khz(125 RBs), and a subband is 6.4Mhz If the subcarrier is
  // 30khz, a RB is 400khz(250 RBs), and a subband is 3.2Mhz
  else if (nof_prb <= 512) {
    return 32;
  } else {
    throw std::runtime_error("Num of physical RBs above 512");
  }
}

static int get_rbg_size(int nof_prb) {
  if (nof_prb <= 10)
    return 1;
  else if (nof_prb <= 26)
    return 2;
  else if (nof_prb <= 63)
    return 3;
  else if (nof_prb <= 110)
    return 4;
  // 5G:
  // https://www.etsi.org/deliver/etsi_ts/138200_138299/138214/15.02.00_60/ts_138214v150200p.pdf
  // : Table 5.1.2.2.1-1 5G:
  // https://www.sharetechnote.com/html/5G/5G_ResourceAllocationType.html
  // (73-144 bandwidth part size => 8 RBs as a RBG) LTE:
  // https://www.sharetechnote.com/html/Handbook_LTE_RAType.html (100 RBs => 4
  // RBs as a RBG)
  else if (nof_prb <= 512)
    return 8;
  else {
    throw std::runtime_error("Num of physical RBs above 512");
  }
}

#endif /* EESM_EFFECTIVE_SINR_H_ */
