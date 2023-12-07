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

#ifndef MULTIPATH_V0_M10_H_
#define MULTIPATH_V0_M10_H_

static float multipath_M10_v_0[3000] = {
    6.18319,     6.18318,    6.18316,    6.18314,    6.18312,    6.18309,
    6.18305,     6.183,      6.18295,    6.1829,     6.18284,    6.18277,
    6.1827,      6.18262,    6.18253,    6.18244,    6.18234,    6.18224,
    6.18213,     6.18202,    6.1819,     6.18177,    6.18164,    6.1815,
    6.18136,     6.18121,    6.18106,    6.18089,    6.18073,    6.18056,
    6.18038,     6.18019,    6.18,       6.17981,    6.1796,     6.1794,
    6.17918,     6.17896,    6.17874,    6.1785,     6.17827,    6.17802,
    6.17777,     6.17752,    6.17726,    6.17699,    6.17672,    6.17644,
    6.17616,     6.17587,    6.17557,    6.17527,    6.17496,    6.17465,
    6.17433,     6.174,      6.17367,    6.17333,    6.17299,    6.17264,
    6.17229,     6.17193,    6.17156,    6.17119,    6.17081,    6.17042,
    6.17003,     6.16964,    6.16924,    6.16883,    6.16841,    6.16799,
    6.16757,     6.16714,    6.1667,     6.16626,    6.16581,    6.16535,
    6.16489,     6.16442,    6.16395,    6.16347,    6.16299,    6.1625,
    6.162,       6.1615,     6.16099,    6.16048,    6.15996,    6.15943,
    6.1589,      6.15836,    6.15781,    6.15727,    6.15671,    6.15615,
    6.15558,     6.15501,    6.15443,    6.15384,    6.15325,    6.15265,
    6.15205,     6.15144,    6.15083,    6.1502,     6.14958,    6.14894,
    6.14831,     6.14766,    6.14701,    6.14635,    6.14569,    6.14502,
    6.14435,     6.14367,    6.14298,    6.14229,    6.14159,    6.14088,
    6.14017,     6.13946,    6.13873,    6.138,      6.13727,    6.13653,
    6.13578,     6.13503,    6.13427,    6.13351,    6.13274,    6.13196,
    6.13118,     6.13039,    6.1296,     6.12879,    6.12799,    6.12718,
    6.12636,     6.12553,    6.1247,     6.12386,    6.12302,    6.12217,
    6.12132,     6.12046,    6.11959,    6.11872,    6.11784,    6.11695,
    6.11606,     6.11517,    6.11426,    6.11335,    6.11244,    6.11152,
    6.11059,     6.10966,    6.10872,    6.10777,    6.10682,    6.10586,
    6.1049,      6.10393,    6.10295,    6.10197,    6.10098,    6.09999,
    6.09899,     6.09798,    6.09697,    6.09595,    6.09492,    6.09389,
    6.09286,     6.09181,    6.09076,    6.08971,    6.08865,    6.08758,
    6.08651,     6.08543,    6.08434,    6.08325,    6.08215,    6.08104,
    6.07993,     6.07882,    6.07769,    6.07656,    6.07543,    6.07429,
    6.07314,     6.07198,    6.07082,    6.06966,    6.06849,    6.06731,
    6.06612,     6.06493,    6.06373,    6.06253,    6.06132,    6.0601,
    6.05888,     6.05765,    6.05642,    6.05518,    6.05393,    6.05268,
    6.05142,     6.05015,    6.04888,    6.0476,     6.04632,    6.04503,
    6.04373,     6.04243,    6.04112,    6.0398,     6.03848,    6.03715,
    6.03581,     6.03447,    6.03313,    6.03177,    6.03041,    6.02905,
    6.02767,     6.02629,    6.02491,    6.02352,    6.02212,    6.02071,
    6.0193,      6.01789,    6.01646,    6.01503,    6.0136,     6.01215,
    6.01071,     6.00925,    6.00779,    6.00632,    6.00485,    6.00337,
    6.00188,     6.00038,    5.99888,    5.99738,    5.99586,    5.99435,
    5.99282,     5.99129,    5.98975,    5.9882,     5.98665,    5.98509,
    5.98353,     5.98196,    5.98038,    5.9788,     5.97721,    5.97561,
    5.97401,     5.9724,     5.97078,    5.96916,    5.96753,    5.96589,
    5.96425,     5.9626,     5.96095,    5.95929,    5.95762,    5.95594,
    5.95426,     5.95257,    5.95088,    5.94918,    5.94747,    5.94576,
    5.94404,     5.94231,    5.94057,    5.93883,    5.93709,    5.93533,
    5.93357,     5.9318,     5.93003,    5.92825,    5.92646,    5.92467,
    5.92287,     5.92106,    5.91925,    5.91743,    5.9156,     5.91377,
    5.91193,     5.91008,    5.90823,    5.90637,    5.9045,     5.90263,
    5.90074,     5.89886,    5.89696,    5.89506,    5.89315,    5.89124,
    5.88932,     5.88739,    5.88546,    5.88352,    5.88157,    5.87961,
    5.87765,     5.87568,    5.87371,    5.87172,    5.86974,    5.86774,
    5.86574,     5.86373,    5.86171,    5.85969,    5.85766,    5.85562,
    5.85358,     5.85153,    5.84947,    5.8474,     5.84533,    5.84325,
    5.84117,     5.83908,    5.83698,    5.83487,    5.83276,    5.83064,
    5.82851,     5.82638,    5.82424,    5.82209,    5.81993,    5.81777,
    5.8156,      5.81343,    5.81124,    5.80905,    5.80686,    5.80465,
    5.80244,     5.80022,    5.798,      5.79576,    5.79352,    5.79128,
    5.78902,     5.78676,    5.7845,     5.78222,    5.77994,    5.77765,
    5.77535,     5.77305,    5.77074,    5.76842,    5.76609,    5.76376,
    5.76142,     5.75908,    5.75672,    5.75436,    5.75199,    5.74962,
    5.74723,     5.74484,    5.74245,    5.74004,    5.73763,    5.73521,
    5.73278,     5.73035,    5.72791,    5.72546,    5.723,      5.72054,
    5.71807,     5.71559,    5.71311,    5.71061,    5.70811,    5.70561,
    5.70309,     5.70057,    5.69804,    5.6955,     5.69296,    5.69041,
    5.68785,     5.68528,    5.6827,     5.68012,    5.67753,    5.67494,
    5.67233,     5.66972,    5.6671,     5.66447,    5.66184,    5.6592,
    5.65655,     5.65389,    5.65122,    5.64855,    5.64587,    5.64318,
    5.64049,     5.63778,    5.63507,    5.63236,    5.62963,    5.6269,
    5.62415,     5.6214,     5.61865,    5.61588,    5.61311,    5.61033,
    5.60754,     5.60475,    5.60194,    5.59913,    5.59631,    5.59349,
    5.59065,     5.58781,    5.58496,    5.5821,     5.57923,    5.57636,
    5.57348,     5.57059,    5.56769,    5.56478,    5.56187,    5.55895,
    5.55602,     5.55308,    5.55014,    5.54718,    5.54422,    5.54125,
    5.53828,     5.53529,    5.5323,     5.52929,    5.52628,    5.52327,
    5.52024,     5.51721,    5.51417,    5.51111,    5.50806,    5.50499,
    5.50192,     5.49883,    5.49574,    5.49264,    5.48953,    5.48642,
    5.48329,     5.48016,    5.47702,    5.47387,    5.47072,    5.46755,
    5.46438,     5.46119,    5.458,      5.4548,     5.4516,     5.44838,
    5.44516,     5.44193,    5.43868,    5.43544,    5.43218,    5.42891,
    5.42564,     5.42235,    5.41906,    5.41576,    5.41245,    5.40914,
    5.40581,     5.40248,    5.39913,    5.39578,    5.39242,    5.38905,
    5.38567,     5.38229,    5.37889,    5.37549,    5.37208,    5.36866,
    5.36523,     5.36179,    5.35834,    5.35489,    5.35142,    5.34795,
    5.34446,     5.34097,    5.33747,    5.33396,    5.33045,    5.32692,
    5.32338,     5.31984,    5.31629,    5.31272,    5.30915,    5.30557,
    5.30198,     5.29838,    5.29478,    5.29116,    5.28753,    5.2839,
    5.28025,     5.2766,     5.27294,    5.26927,    5.26559,    5.2619,
    5.2582,      5.25449,    5.25077,    5.24705,    5.24331,    5.23957,
    5.23581,     5.23205,    5.22828,    5.22449,    5.2207,     5.2169,
    5.21309,     5.20927,    5.20544,    5.2016,     5.19775,    5.1939,
    5.19003,     5.18615,    5.18227,    5.17837,    5.17446,    5.17055,
    5.16662,     5.16269,    5.15875,    5.15479,    5.15083,    5.14685,
    5.14287,     5.13888,    5.13488,    5.13086,    5.12684,    5.12281,
    5.11877,     5.11472,    5.11066,    5.10659,    5.1025,     5.09841,
    5.09431,     5.0902,     5.08608,    5.08195,    5.07781,    5.07366,
    5.06949,     5.06532,    5.06114,    5.05695,    5.05275,    5.04854,
    5.04431,     5.04008,    5.03584,    5.03158,    5.02732,    5.02305,
    5.01876,     5.01447,    5.01017,    5.00585,    5.00152,    4.99719,
    4.99284,     4.98849,    4.98412,    4.97974,    4.97535,    4.97095,
    4.96654,     4.96212,    4.95769,    4.95325,    4.9488,     4.94433,
    4.93986,     4.93538,    4.93088,    4.92637,    4.92186,    4.91733,
    4.91279,     4.90824,    4.90368,    4.89911,    4.89453,    4.88993,
    4.88533,     4.88071,    4.87609,    4.87145,    4.8668,     4.86214,
    4.85747,     4.85279,    4.84809,    4.84339,    4.83867,    4.83395,
    4.82921,     4.82446,    4.8197,     4.81492,    4.81014,    4.80534,
    4.80054,     4.79572,    4.79089,    4.78605,    4.78119,    4.77633,
    4.77145,     4.76656,    4.76166,    4.75675,    4.75183,    4.74689,
    4.74194,     4.73698,    4.73201,    4.72703,    4.72204,    4.71703,
    4.71201,     4.70698,    4.70194,    4.69688,    4.69182,    4.68674,
    4.68165,     4.67654,    4.67143,    4.6663,     4.66116,    4.65601,
    4.65084,     4.64566,    4.64048,    4.63527,    4.63006,    4.62483,
    4.61959,     4.61434,    4.60908,    4.6038,     4.59851,    4.59321,
    4.58789,     4.58256,    4.57722,    4.57187,    4.5665,     4.56112,
    4.55573,     4.55032,    4.5449,     4.53947,    4.53403,    4.52857,
    4.5231,      4.51762,    4.51212,    4.50661,    4.50109,    4.49555,
    4.49,        4.48444,    4.47886,    4.47327,    4.46766,    4.46205,
    4.45642,     4.45077,    4.44511,    4.43944,    4.43375,    4.42806,
    4.42234,     4.41662,    4.41087,    4.40512,    4.39935,    4.39357,
    4.38777,     4.38196,    4.37614,    4.3703,     4.36444,    4.35858,
    4.3527,      4.3468,     4.34089,    4.33496,    4.32903,    4.32307,
    4.3171,      4.31112,    4.30512,    4.29911,    4.29309,    4.28705,
    4.28099,     4.27492,    4.26883,    4.26273,    4.25662,    4.25049,
    4.24434,     4.23818,    4.23201,    4.22582,    4.21961,    4.21339,
    4.20715,     4.2009,     4.19464,    4.18835,    4.18206,    4.17574,
    4.16941,     4.16307,    4.15671,    4.15033,    4.14394,    4.13753,
    4.13111,     4.12467,    4.11822,    4.11174,    4.10526,    4.09875,
    4.09223,     4.0857,     4.07915,    4.07258,    4.06599,    4.05939,
    4.05278,     4.04614,    4.03949,    4.03282,    4.02614,    4.01944,
    4.01272,     4.00599,    3.99924,    3.99247,    3.98568,    3.97888,
    3.97206,     3.96523,    3.95837,    3.9515,     3.94461,    3.93771,
    3.93079,     3.92385,    3.91689,    3.90991,    3.90292,    3.89591,
    3.88888,     3.88183,    3.87477,    3.86769,    3.86058,    3.85347,
    3.84633,     3.83917,    3.832,      3.82481,    3.8176,     3.81037,
    3.80312,     3.79586,    3.78857,    3.78127,    3.77395,    3.76661,
    3.75925,     3.75187,    3.74447,    3.73705,    3.72962,    3.72216,
    3.71468,     3.70719,    3.69968,    3.69214,    3.68459,    3.67702,
    3.66942,     3.66181,    3.65418,    3.64653,    3.63886,    3.63116,
    3.62345,     3.61572,    3.60796,    3.60019,    3.5924,     3.58458,
    3.57675,     3.56889,    3.56101,    3.55311,    3.5452,     3.53726,
    3.52929,     3.52131,    3.51331,    3.50528,    3.49724,    3.48917,
    3.48108,     3.47297,    3.46483,    3.45668,    3.4485,     3.4403,
    3.43208,     3.42384,    3.41557,    3.40728,    3.39897,    3.39064,
    3.38229,     3.37391,    3.36551,    3.35708,    3.34863,    3.34016,
    3.33167,     3.32315,    3.31461,    3.30605,    3.29746,    3.28885,
    3.28021,     3.27155,    3.26287,    3.25416,    3.24543,    3.23668,
    3.2279,      3.21909,    3.21026,    3.20141,    3.19253,    3.18363,
    3.1747,      3.16574,    3.15677,    3.14776,    3.13873,    3.12968,
    3.12059,     3.11149,    3.10236,    3.0932,     3.08401,    3.0748,
    3.06556,     3.0563,     3.04701,    3.03769,    3.02835,    3.01898,
    3.00958,     3.00015,    2.9907,     2.98122,    2.97171,    2.96218,
    2.95262,     2.94302,    2.93341,    2.92376,    2.91408,    2.90438,
    2.89465,     2.88488,    2.87509,    2.86528,    2.85543,    2.84555,
    2.83564,     2.82571,    2.81574,    2.80575,    2.79572,    2.78566,
    2.77558,     2.76546,    2.75531,    2.74514,    2.73493,    2.72469,
    2.71442,     2.70412,    2.69378,    2.68342,    2.67302,    2.66259,
    2.65213,     2.64164,    2.63111,    2.62056,    2.60996,    2.59934,
    2.58869,     2.578,      2.56727,    2.55652,    2.54573,    2.5349,
    2.52404,     2.51315,    2.50222,    2.49126,    2.48027,    2.46924,
    2.45817,     2.44707,    2.43593,    2.42476,    2.41355,    2.4023,
    2.39102,     2.37971,    2.36835,    2.35696,    2.34553,    2.33407,
    2.32256,     2.31102,    2.29944,    2.28783,    2.27617,    2.26448,
    2.25275,     2.24098,    2.22917,    2.21732,    2.20543,    2.1935,
    2.18153,     2.16952,    2.15747,    2.14538,    2.13325,    2.12108,
    2.10886,     2.09661,    2.08431,    2.07197,    2.05959,    2.04716,
    2.0347,      2.02219,    2.00963,    1.99703,    1.98439,    1.97171,
    1.95898,     1.9462,     1.93338,    1.92051,    1.9076,     1.89465,
    1.88164,     1.86859,    1.8555,     1.84235,    1.82916,    1.81593,
    1.80264,     1.78931,    1.77592,    1.76249,    1.74901,    1.73548,
    1.7219,      1.70827,    1.69459,    1.68086,    1.66708,    1.65325,
    1.63936,     1.62543,    1.61144,    1.5974,     1.5833,     1.56915,
    1.55495,     1.54069,    1.52638,    1.51202,    1.4976,     1.48312,
    1.46859,     1.454,      1.43936,    1.42465,    1.40989,    1.39507,
    1.3802,      1.36526,    1.35027,    1.33521,    1.3201,     1.30493,
    1.28969,     1.27439,    1.25903,    1.24361,    1.22813,    1.21258,
    1.19697,     1.1813,     1.16556,    1.14976,    1.13389,    1.11795,
    1.10195,     1.08588,    1.06974,    1.05354,    1.03727,    1.02093,
    1.00451,     0.988032,   0.971479,   0.954855,   0.93816,    0.921392,
    0.904551,    0.887637,   0.870648,   0.853585,   0.836447,   0.819232,
    0.801941,    0.784573,   0.767128,   0.749603,   0.732,      0.714317,
    0.696554,    0.678709,   0.660783,   0.642774,   0.624683,   0.606507,
    0.588247,    0.569901,   0.55147,    0.532952,   0.514346,   0.495652,
    0.476869,    0.457996,   0.439032,   0.419977,   0.40083,    0.381589,
    0.362255,    0.342825,   0.3233,     0.303679,   0.28396,    0.264142,
    0.244226,    0.224209,   0.204091,   0.183871,   0.163548,   0.14312,
    0.122588,    0.101949,   0.0812038,  0.0603501,  0.0393872,  0.0183139,
    -0.00287071, -0.0241679, -0.0455786, -0.0671042, -0.0887457, -0.110504,
    -0.132381,   -0.154378,  -0.176495,  -0.198734,  -0.221097,  -0.243584,
    -0.266197,   -0.288937,  -0.311806,  -0.334805,  -0.357935,  -0.381198,
    -0.404596,   -0.428129,  -0.451799,  -0.475609,  -0.499558,  -0.52365,
    -0.547884,   -0.572265,  -0.596792,  -0.621467,  -0.646293,  -0.67127,
    -0.696402,   -0.721689,  -0.747133,  -0.772737,  -0.798502,  -0.82443,
    -0.850524,   -0.876784,  -0.903214,  -0.929815,  -0.956589,  -0.98354,
    -1.01067,    -1.03798,   -1.06547,   -1.09314,   -1.12101,   -1.14906,
    -1.1773,     -1.20574,   -1.23438,   -1.26322,   -1.29225,   -1.3215,
    -1.35095,    -1.38062,   -1.41049,   -1.44059,   -1.4709,    -1.50144,
    -1.53221,    -1.5632,    -1.59443,   -1.6259,    -1.6576,    -1.68955,
    -1.72174,    -1.75419,   -1.78689,   -1.81985,   -1.85307,   -1.88656,
    -1.92032,    -1.95435,   -1.98867,   -2.02327,   -2.05816,   -2.09334,
    -2.12882,    -2.1646,    -2.20069,   -2.23709,   -2.27382,   -2.31086,
    -2.34824,    -2.38595,   -2.42401,   -2.46241,   -2.50117,   -2.54029,
    -2.57977,    -2.61963,   -2.65987,   -2.7005,    -2.74152,   -2.78295,
    -2.82479,    -2.86705,   -2.90973,   -2.95285,   -2.99641,   -3.04043,
    -3.08491,    -3.12987,   -3.17531,   -3.22124,   -3.26767,   -3.31462,
    -3.3621,     -3.41011,   -3.45867,   -3.5078,    -3.5575,    -3.60779,
    -3.65868,    -3.71019,   -3.76233,   -3.81512,   -3.86857,   -3.9227,
    -3.97753,    -4.03308,   -4.08936,   -4.14639,   -4.20419,   -4.26279,
    -4.32221,    -4.38247,   -4.44359,   -4.50559,   -4.56851,   -4.63238,
    -4.69721,    -4.76304,   -4.82989,   -4.89781,   -4.96683,   -5.03698,
    -5.10829,    -5.18081,   -5.25459,   -5.32965,   -5.40605,   -5.48384,
    -5.56307,    -5.64378,   -5.72605,   -5.80992,   -5.89546,   -5.98274,
    -6.07183,    -6.16281,   -6.25575,   -6.35075,   -6.44789,   -6.54728,
    -6.64902,    -6.75322,   -6.86001,   -6.96951,   -7.08186,   -7.19723,
    -7.31577,    -7.43766,   -7.56309,   -7.69229,   -7.82547,   -7.96289,
    -8.10484,    -8.25161,   -8.40354,   -8.56102,   -8.72445,   -8.89431,
    -9.07112,    -9.25547,   -9.44803,   -9.64956,   -9.86094,   -10.0832,
    -10.3175,    -10.5652,   -10.8279,   -11.1076,   -11.4066,   -11.7278,
    -12.0747,    -12.4518,   -12.8649,   -13.3216,   -13.832,    -14.4107,
    -15.0786,    -15.8684,   -16.8346,   -18.0795,   -19.8317,   -22.8168,
    -42.2375,    -22.9164,   -19.8805,   -18.1108,   -16.8569,   -15.8849,
    -15.0911,    -14.4201,   -13.8389,   -13.3263,   -12.8679,   -12.4532,
    -12.0745,    -11.7263,   -11.4038,   -11.1036,   -10.8228,   -10.559,
    -10.3103,    -10.075,    -9.85176,   -9.63943,   -9.43697,   -9.2435,
    -9.05826,    -8.88057,   -8.70985,   -8.54556,   -8.38724,   -8.23447,
    -8.08688,    -7.94412,   -7.80588,   -7.67189,   -7.5419,    -7.41567,
    -7.29299,    -7.17366,   -7.05752,   -6.94438,   -6.83411,   -6.72655,
    -6.62158,    -6.51908,   -6.41893,   -6.32102,   -6.22526,   -6.13156,
    -6.03983,    -5.94999,   -5.86195,   -5.77566,   -5.69104,   -5.60803,
    -5.52657,    -5.4466,    -5.36807,   -5.29092,   -5.21512,   -5.1406,
    -5.06734,    -4.99529,   -4.9244,    -4.85465,   -4.786,     -4.7184,
    -4.65184,    -4.58627,   -4.52168,   -4.45803,   -4.39529,   -4.33344,
    -4.27245,    -4.2123,    -4.15297,   -4.09443,   -4.03667,   -3.97967,
    -3.92339,    -3.86784,   -3.81298,   -3.7588,    -3.70528,   -3.65241,
    -3.60018,    -3.54856,   -3.49755,   -3.44712,   -3.39727,   -3.34798,
    -3.29924,    -3.25105,   -3.20337,   -3.15621,   -3.10956,   -3.0634,
    -3.01772,    -2.97252,   -2.92777,   -2.88349,   -2.83964,   -2.79624,
    -2.75326,    -2.7107,    -2.66855,   -2.6268,    -2.58545,   -2.54449,
    -2.50391,    -2.46371,   -2.42387,   -2.38439,   -2.34526,   -2.30649,
    -2.26805,    -2.22995,   -2.19219,   -2.15474,   -2.11762,   -2.08081,
    -2.0443,     -2.0081,    -1.9722,    -1.93659,   -1.90127,   -1.86624,
    -1.83148,    -1.797,     -1.7628,    -1.72885,   -1.69517,   -1.66175,
    -1.62859,    -1.59567,   -1.563,     -1.53058,   -1.4984,    -1.46645,
    -1.43473,    -1.40325,   -1.37199,   -1.34095,   -1.31014,   -1.27954,
    -1.24916,    -1.21898,   -1.18902,   -1.15926,   -1.1297,    -1.10035,
    -1.07119,    -1.04222,   -1.01345,   -0.984868,  -0.956473,  -0.928262,
    -0.900234,   -0.872386,  -0.844715,  -0.817221,  -0.789899,  -0.762749,
    -0.735769,   -0.708955,  -0.682307,  -0.655822,  -0.629497,  -0.603333,
    -0.577325,   -0.551473,  -0.525775,  -0.500229,  -0.474833,  -0.449585,
    -0.424485,   -0.399529,  -0.374717,  -0.350047,  -0.325517,  -0.301125,
    -0.276871,   -0.252753,  -0.228769,  -0.204917,  -0.181197,  -0.157607,
    -0.134146,   -0.110811,  -0.0876027, -0.0645187, -0.0415579, -0.0187192,
    0.00399889,  0.0265975,  0.0490778,  0.0714411,  0.0936886,  0.115821,
    0.13784,     0.159747,   0.181543,   0.203228,   0.224804,   0.246272,
    0.267633,    0.288888,   0.310037,   0.331083,   0.352026,   0.372867,
    0.393607,    0.414247,   0.434787,   0.455229,   0.475574,   0.495822,
    0.515975,    0.536033,   0.555997,   0.575868,   0.595647,   0.615334,
    0.634931,    0.654438,   0.673857,   0.693187,   0.712429,   0.731585,
    0.750655,    0.76964,    0.78854,    0.807356,   0.82609,    0.844741,
    0.86331,     0.881799,   0.900207,   0.918536,   0.936785,   0.954957,
    0.97305,     0.991067,   1.00901,    1.02687,    1.04466,    1.06238,
    1.08002,     1.09758,    1.11508,    1.1325,     1.14985,    1.16713,
    1.18434,     1.20148,    1.21855,    1.23555,    1.25248,    1.26934,
    1.28614,     1.30287,    1.31953,    1.33613,    1.35266,    1.36913,
    1.38553,     1.40187,    1.41814,    1.43436,    1.45051,    1.46659,
    1.48262,     1.49858,    1.51449,    1.53033,    1.54611,    1.56184,
    1.5775,      1.59311,    1.60866,    1.62415,    1.63958,    1.65496,
    1.67028,     1.68554,    1.70075,    1.7159,     1.73099,    1.74604,
    1.76102,     1.77596,    1.79084,    1.80566,    1.82044,    1.83516,
    1.84983,     1.86444,    1.87901,    1.89352,    1.90799,    1.9224,
    1.93676,     1.95107,    1.96534,    1.97955,    1.99371,    2.00783,
    2.0219,      2.03592,    2.04989,    2.06381,    2.07769,    2.09152,
    2.10531,     2.11904,    2.13274,    2.14638,    2.15998,    2.17354,
    2.18705,     2.20052,    2.21394,    2.22732,    2.24065,    2.25394,
    2.26719,     2.28039,    2.29355,    2.30667,    2.31975,    2.33279,
    2.34578,     2.35873,    2.37164,    2.38451,    2.39734,    2.41013,
    2.42287,     2.43558,    2.44825,    2.46088,    2.47346,    2.48601,
    2.49852,     2.51099,    2.52343,    2.53582,    2.54818,    2.5605,
    2.57278,     2.58502,    2.59723,    2.6094,     2.62153,    2.63362,
    2.64568,     2.65771,    2.66969,    2.68164,    2.69356,    2.70544,
    2.71728,     2.72909,    2.74086,    2.7526,     2.76431,    2.77598,
    2.78761,     2.79921,    2.81078,    2.82232,    2.83382,    2.84529,
    2.85672,     2.86812,    2.87949,    2.89082,    2.90213,    2.9134,
    2.92464,     2.93584,    2.94702,    2.95816,    2.96927,    2.98035,
    2.9914,      3.00242,    3.0134,     3.02436,    3.03529,    3.04618,
    3.05704,     3.06788,    3.07868,    3.08946,    3.1002,     3.11092,
    3.1216,      3.13226,    3.14289,    3.15348,    3.16405,    3.17459,
    3.1851,      3.19559,    3.20604,    3.21647,    3.22687,    3.23724,
    3.24758,     3.25789,    3.26818,    3.27844,    3.28867,    3.29888,
    3.30906,     3.31921,    3.32933,    3.33943,    3.3495,     3.35954,
    3.36956,     3.37955,    3.38952,    3.39946,    3.40937,    3.41926,
    3.42912,     3.43896,    3.44877,    3.45855,    3.46831,    3.47805,
    3.48776,     3.49745,    3.50711,    3.51674,    3.52635,    3.53594,
    3.5455,      3.55504,    3.56455,    3.57404,    3.58351,    3.59295,
    3.60237,     3.61177,    3.62114,    3.63048,    3.63981,    3.64911,
    3.65839,     3.66764,    3.67687,    3.68608,    3.69527,    3.70443,
    3.71357,     3.72269,    3.73179,    3.74086,    3.74991,    3.75894,
    3.76795,     3.77693,    3.78589,    3.79483,    3.80375,    3.81265,
    3.82153,     3.83038,    3.83922,    3.84803,    3.85682,    3.86559,
    3.87434,     3.88306,    3.89177,    3.90046,    3.90912,    3.91777,
    3.92639,     3.93499,    3.94358,    3.95214,    3.96068,    3.9692,
    3.97771,     3.98619,    3.99465,    4.00309,    4.01152,    4.01992,
    4.0283,      4.03666,    4.04501,    4.05333,    4.06164,    4.06992,
    4.07819,     4.08644,    4.09467,    4.10288,    4.11107,    4.11924,
    4.12739,     4.13553,    4.14364,    4.15174,    4.15982,    4.16788,
    4.17592,     4.18394,    4.19194,    4.19993,    4.2079,     4.21585,
    4.22378,     4.2317,     4.23959,    4.24747,    4.25533,    4.26317,
    4.271,       4.27881,    4.2866,     4.29437,    4.30213,    4.30986,
    4.31759,     4.32529,    4.33298,    4.34064,    4.3483,     4.35593,
    4.36355,     4.37115,    4.37874,    4.3863,     4.39385,    4.40139,
    4.40891,     4.41641,    4.42389,    4.43136,    4.43881,    4.44625,
    4.45367,     4.46107,    4.46846,    4.47583,    4.48318,    4.49052,
    4.49785,     4.50515,    4.51244,    4.51972,    4.52698,    4.53422,
    4.54145,     4.54866,    4.55586,    4.56304,    4.57021,    4.57736,
    4.58449,     4.59161,    4.59872,    4.60581,    4.61288,    4.61994,
    4.62699,     4.63402,    4.64103,    4.64803,    4.65501,    4.66198,
    4.66894,     4.67588,    4.6828,     4.68971,    4.69661,    4.70349,
    4.71036,     4.71721,    4.72405,    4.73087,    4.73768,    4.74448,
    4.75126,     4.75802,    4.76478,    4.77151,    4.77824,    4.78495,
    4.79164,     4.79833,    4.805,      4.81165,    4.81829,    4.82492,
    4.83153,     4.83813,    4.84471,    4.85129,    4.85784,    4.86439,
    4.87092,     4.87744,    4.88394,    4.89043,    4.89691,    4.90337,
    4.90982,     4.91626,    4.92269,    4.9291,     4.93549,    4.94188,
    4.94825,     4.95461,    4.96096,    4.96729,    4.97361,    4.97992,
    4.98621,     4.99249,    4.99876,    5.00501,    5.01126,    5.01749,
    5.02371,     5.02991,    5.0361,     5.04228,    5.04845,    5.05461,
    5.06075,     5.06688,    5.073,      5.0791,     5.08519,    5.09128,
    5.09734,     5.1034,     5.10944,    5.11548,    5.1215,     5.1275,
    5.1335,      5.13948,    5.14546,    5.15142,    5.15736,    5.1633,
    5.16922,     5.17514,    5.18104,    5.18693,    5.1928,     5.19867,
    5.20452,     5.21036,    5.2162,     5.22201,    5.22782,    5.23362,
    5.2394,      5.24517,    5.25094,    5.25669,    5.26242,    5.26815,
    5.27387,     5.27957,    5.28527,    5.29095,    5.29662,    5.30228,
    5.30793,     5.31357,    5.31919,    5.32481,    5.33041,    5.336,
    5.34159,     5.34716,    5.35272,    5.35827,    5.36381,    5.36933,
    5.37485,     5.38036,    5.38585,    5.39134,    5.39681,    5.40227,
    5.40773,     5.41317,    5.4186,     5.42402,    5.42943,    5.43483,
    5.44022,     5.4456,     5.45097,    5.45632,    5.46167,    5.46701,
    5.47233,     5.47765,    5.48295,    5.48825,    5.49354,    5.49881,
    5.50407,     5.50933,    5.51457,    5.51981,    5.52503,    5.53025,
    5.53545,     5.54064,    5.54583,    5.551,      5.55616,    5.56132,
    5.56646,     5.57159,    5.57672,    5.58183,    5.58694,    5.59203,
    5.59712,     5.60219,    5.60725,    5.61231,    5.61736,    5.62239,
    5.62742,     5.63243,    5.63744,    5.64244,    5.64742,    5.6524,
    5.65737,     5.66233,    5.66728,    5.67222,    5.67715,    5.68207,
    5.68698,     5.69188,    5.69677,    5.70166,    5.70653,    5.7114,
    5.71625,     5.7211,     5.72593,    5.73076,    5.73558,    5.74039,
    5.74519,     5.74998,    5.75476,    5.75953,    5.76429,    5.76905,
    5.77379,     5.77853,    5.78326,    5.78797,    5.79268,    5.79738,
    5.80207,     5.80675,    5.81143,    5.81609,    5.82075,    5.82539,
    5.83003,     5.83466,    5.83928,    5.84389,    5.84849,    5.85308,
    5.85767,     5.86224,    5.86681,    5.87137,    5.87592,    5.88046,
    5.88499,     5.88952,    5.89403,    5.89854,    5.90303,    5.90752,
    5.912,       5.91648,    5.92094,    5.9254,     5.92984,    5.93428,
    5.93871,     5.94313,    5.94754,    5.95195,    5.95635,    5.96073,
    5.96511,     5.96948,    5.97385,    5.9782,     5.98255,    5.98688,
    5.99121,     5.99554,    5.99985,    6.00415,    6.00845,    6.01274,
    6.01702,     6.02129,    6.02556,    6.02981,    6.03406,    6.0383,
    6.04253,     6.04676,    6.05097,    6.05518,    6.05938,    6.06357,
    6.06775,     6.07193,    6.0761,     6.08026,    6.08441,    6.08855,
    6.09269,     6.09682,    6.10094,    6.10505,    6.10915,    6.11325,
    6.11734,     6.12142,    6.12549,    6.12956,    6.13362,    6.13767,
    6.14171,     6.14574,    6.14977,    6.15379,    6.1578,     6.16181,
    6.1658,      6.16979,    6.17377,    6.17774,    6.18171,    6.18567,
    6.18962,     6.19356,    6.1975,     6.20143,    6.20535,    6.20926,
    6.21316,     6.21706,    6.22095,    6.22484,    6.22871,    6.23258,
    6.23644,     6.24029,    6.24414,    6.24798,    6.25181,    6.25563,
    6.25945,     6.26326,    6.26706,    6.27086,    6.27464,    6.27842,
    6.2822,      6.28596,    6.28972,    6.29347,    6.29722,    6.30095,
    6.30468,     6.3084,     6.31212,    6.31583,    6.31953,    6.32322,
    6.32691,     6.33059,    6.33426,    6.33792,    6.34158,    6.34523,
    6.34888,     6.35251,    6.35614,    6.35977,    6.36338,    6.36699,
    6.37059,     6.37419,    6.37778,    6.38136,    6.38493,    6.3885,
    6.39206,     6.39561,    6.39916,    6.4027,     6.40623,    6.40976,
    6.41328,     6.41679,    6.42029,    6.42379,    6.42728,    6.43077,
    6.43424,     6.43772,    6.44118,    6.44464,    6.44809,    6.45153,
    6.45497,     6.4584,     6.46182,    6.46524,    6.46865,    6.47206,
    6.47545,     6.47884,    6.48223,    6.4856,     6.48897,    6.49234,
    6.4957,      6.49905,    6.50239,    6.50573,    6.50906,    6.51238,
    6.5157,      6.51901,    6.52232,    6.52561,    6.52891,    6.53219,
    6.53547,     6.53874,    6.54201,    6.54527,    6.54852,    6.55177,
    6.55501,     6.55824,    6.56147,    6.56469,    6.5679,     6.57111,
    6.57431,     6.5775,     6.58069,    6.58387,    6.58705,    6.59022,
    6.59338,     6.59654,    6.59969,    6.60283,    6.60597,    6.6091,
    6.61223,     6.61535,    6.61846,    6.62157,    6.62467,    6.62776,
    6.63085,     6.63393,    6.63701,    6.64007,    6.64314,    6.64619,
    6.64924,     6.65229,    6.65533,    6.65836,    6.66139,    6.6644,
    6.66742,     6.67043,    6.67343,    6.67642,    6.67941,    6.68239,
    6.68537,     6.68834,    6.69131,    6.69427,    6.69722,    6.70016,
    6.7031,      6.70604,    6.70897,    6.71189,    6.71481,    6.71772,
    6.72062,     6.72352,    6.72641,    6.7293,     6.73218,    6.73505,
    6.73792,     6.74079,    6.74364,    6.74649,    6.74934,    6.75218,
    6.75501,     6.75784,    6.76066,    6.76347,    6.76628,    6.76909,
    6.77189,     6.77468,    6.77746,    6.78024,    6.78302,    6.78579,
    6.78855,     6.79131,    6.79406,    6.7968,     6.79954,    6.80228,
    6.805,       6.80773,    6.81044,    6.81315,    6.81586,    6.81856,
    6.82125,     6.82394,    6.82662,    6.8293,     6.83197,    6.83463,
    6.83729,     6.83995,    6.84259,    6.84524,    6.84787,    6.8505,
    6.85313,     6.85575,    6.85836,    6.86097,    6.86357,    6.86617,
    6.86876,     6.87135,    6.87393,    6.8765,     6.87907,    6.88164,
    6.88419,     6.88675,    6.88929,    6.89183,    6.89437,    6.8969,
    6.89942,     6.90194,    6.90446,    6.90697,    6.90947,    6.91196,
    6.91446,     6.91694,    6.91942,    6.9219,     6.92437,    6.92683,
    6.92929,     6.93174,    6.93419,    6.93663,    6.93907,    6.9415,
    6.94393,     6.94635,    6.94876,    6.95117,    6.95358,    6.95597,
    6.95837,     6.96076,    6.96314,    6.96552,    6.96789,    6.97025,
    6.97262,     6.97497,    6.97732,    6.97967,    6.98201,    6.98434,
    6.98667,     6.98899,    6.99131,    6.99363,    6.99593,    6.99824,
    7.00053,     7.00282,    7.00511,    7.00739,    7.00967,    7.01194,
    7.01421,     7.01647,    7.01872,    7.02097,    7.02322,    7.02545,
    7.02769,     7.02992,    7.03214,    7.03436,    7.03657,    7.03878,
    7.04098,     7.04318,    7.04537,    7.04756,    7.04974,    7.05192,
    7.05409,     7.05626,    7.05842,    7.06057,    7.06273,    7.06487,
    7.06701,     7.06915,    7.07128,    7.0734,     7.07552,    7.07764,
    7.07975,     7.08185,    7.08395,    7.08605,    7.08814,    7.09022,
    7.0923,      7.09438,    7.09645,    7.09851,    7.10057,    7.10262,
    7.10467,     7.10672,    7.10875,    7.11079,    7.11282,    7.11484,
    7.11686,     7.11887,    7.12088,    7.12289,    7.12488,    7.12688,
    7.12887,     7.13085,    7.13283,    7.1348,     7.13677,    7.13873,
    7.14069,     7.14265,    7.1446,     7.14654,    7.14848,    7.15041,
    7.15234,     7.15426,    7.15618,    7.1581,     7.16001,    7.16191,
    7.16381,     7.16571,    7.16759,    7.16948,    7.17136,    7.17323,
    7.1751,      7.17697,    7.17883,    7.18068,    7.18253,    7.18438,
    7.18622,     7.18806,    7.18989,    7.19171,    7.19353,    7.19535,
    7.19716,     7.19897,    7.20077,    7.20257,    7.20436,    7.20614,
    7.20793,     7.2097,     7.21148,    7.21324,    7.21501,    7.21677,
    7.21852,     7.22027,    7.22201,    7.22375,    7.22549,    7.22721,
    7.22894,     7.23066,    7.23237,    7.23408,    7.23579,    7.23749,
    7.23919,     7.24088,    7.24256,    7.24425,    7.24592,    7.2476,
    7.24926,     7.25093,    7.25258,    7.25424,    7.25589,    7.25753,
    7.25917,     7.2608,     7.26243,    7.26406,    7.26568,    7.26729,
    7.26891,     7.27051,    7.27211,    7.27371,    7.2753,     7.27689,
    7.27847,     7.28005,    7.28162,    7.28319,    7.28476,    7.28632,
    7.28787,     7.28942,    7.29097,    7.29251,    7.29404,    7.29558,
    7.2971,      7.29862,    7.30014,    7.30166,    7.30316,    7.30467,
    7.30617,     7.30766,    7.30915,    7.31064,    7.31212,    7.31359,
    7.31506,     7.31653,    7.31799,    7.31945,    7.32091,    7.32235,
    7.3238,      7.32524,    7.32667,    7.3281,     7.32953,    7.33095,
    7.33237,     7.33378,    7.33518,    7.33659,    7.33799,    7.33938,
    7.34077,     7.34215,    7.34353,    7.34491,    7.34628,    7.34765,
    7.34901,     7.35037,    7.35172,    7.35307,    7.35441,    7.35575,
    7.35709,     7.35842,    7.35974,    7.36106,    7.36238,    7.36369,
    7.365,       7.3663,     7.3676,     7.3689,     7.37019,    7.37147,
    7.37275,     7.37403,    7.3753,     7.37657,    7.37783,    7.37909,
    7.38035,     7.38159,    7.38284,    7.38408,    7.38532,    7.38655,
    7.38778,     7.389,      7.39022,    7.39143,    7.39264,    7.39385,
    7.39505,     7.39624,    7.39744,    7.39862,    7.39981,    7.40098,
    7.40216,     7.40333,    7.40449,    7.40565,    7.40681,    7.40796,
    7.40911,     7.41025,    7.41139,    7.41253,    7.41366,    7.41478,
    7.41591,     7.41702,    7.41813,    7.41924,    7.42035,    7.42145,
    7.42254,     7.42363,    7.42472,    7.4258,     7.42688,    7.42795,
    7.42902,     7.43008,    7.43114,    7.4322,     7.43325,    7.4343,
    7.43534,     7.43638,    7.43741,    7.43844,    7.43947,    7.44049,
    7.44151,     7.44252,    7.44353,    7.44453,    7.44553,    7.44652,
    7.44751,     7.4485,     7.44948,    7.45046,    7.45143,    7.4524,
    7.45337,     7.45433,    7.45528,    7.45624,    7.45718,    7.45813,
    7.45907,     7.46,       7.46093,    7.46186,    7.46278,    7.4637,
    7.46461,     7.46552,    7.46642,    7.46732,    7.46822,    7.46911,
    7.47,        7.47088,    7.47176,    7.47263,    7.4735,     7.47437,
    7.47523,     7.47609,    7.47694,    7.47779,    7.47864,    7.47948,
    7.48031,     7.48115,    7.48197,    7.4828,     7.48362,    7.48443,
    7.48524,     7.48605,    7.48685,    7.48765,    7.48844,    7.48923,
    7.49002,     7.4908,     7.49157,    7.49235,    7.49312,    7.49388,
    7.49464,     7.49539,    7.49615,    7.49689,    7.49764,    7.49838,
    7.49911,     7.49984,    7.50057,    7.50129,    7.50201,    7.50272,
    7.50343,     7.50413,    7.50483,    7.50553,    7.50622,    7.50691,
    7.5076,      7.50828,    7.50895,    7.50962,    7.51029,    7.51095,
    7.51161,     7.51227,    7.51292,    7.51356,    7.51421,    7.51484,
    7.51548,     7.51611,    7.51673,    7.51735,    7.51797,    7.51858,
    7.51919,     7.5198,     7.5204,     7.52099,    7.52158,    7.52217,
    7.52276,     7.52334,    7.52391,    7.52448,    7.52505,    7.52561,
    7.52617,     7.52672,    7.52728,    7.52782,    7.52836,    7.5289,
    7.52944,     7.52996,    7.53049,    7.53101,    7.53153,    7.53204,
    7.53255,     7.53305,    7.53356,    7.53405,    7.53454,    7.53503,
    7.53552,     7.536,      7.53647,    7.53694,    7.53741,    7.53787,
    7.53833,     7.53879,    7.53924,    7.53969,    7.54013,    7.54057,
    7.541,       7.54143,    7.54186,    7.54228,    7.5427,     7.54311,
    7.54352,     7.54393,    7.54433,    7.54473,    7.54512,    7.54551,
    7.54589,     7.54627,    7.54665,    7.54702,    7.54739,    7.54776,
    7.54812,     7.54847,    7.54882,    7.54917,    7.54952,    7.54986,
    7.55019,     7.55052,    7.55085,    7.55117,    7.55149,    7.55181,
    7.55212,     7.55243,    7.55273,    7.55303,    7.55332,    7.55361,
    7.5539,      7.55418,    7.55446,    7.55474,    7.55501,    7.55527,
    7.55553,     7.55579,    7.55604,    7.55629,    7.55654,    7.55678,
    7.55702,     7.55725,    7.55748,    7.55771,    7.55793,    7.55814,
    7.55836,     7.55857,    7.55877,    7.55897,    7.55917,    7.55936,
    7.55955,     7.55973,    7.55991,    7.56009,    7.56026,    7.56043,
    7.56059,     7.56075,    7.56091,    7.56106,    7.56121,    7.56135,
    7.56149,     7.56163,    7.56176,    7.56189,    7.56201,    7.56213,
    7.56224,     7.56236,    7.56246,    7.56257,    7.56266,    7.56276,
    7.56285,     7.56294,    7.56302,    7.5631,     7.56317,    7.56324,
    7.56331,     7.56337,    7.56343,    7.56348,    7.56353,    7.56358,
    7.56362,     7.56366,    7.56369,    7.56372,    7.56375,    7.56377,
    7.56379,     7.5638,     7.56381,    7.56382,    7.56382,    7.56381,
    7.56381,     7.5638,     7.56378,    7.56376,    7.56374,    7.56371,
    7.56368,     7.56365,    7.56361,    7.56356,    7.56352,    7.56347,
    7.56341,     7.56335,    7.56329,    7.56322,    7.56315,    7.56307,
    7.56299,     7.56291,    7.56282,    7.56273,    7.56263,    7.56253,
    7.56243,     7.56232,    7.56221,    7.56209,    7.56197,    7.56185};

#endif /* MULTIPATH_V0_M10_H_ */
