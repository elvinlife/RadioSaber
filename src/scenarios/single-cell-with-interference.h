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

#include "../channel/LteChannel.h"
#include "../phy/enb-lte-phy.h"
#include "../phy/ue-lte-phy.h"
#include "../core/spectrum/bandwidth-manager.h"
#include "../networkTopology/Cell.h"
#include "../protocolStack/packet/packet-burst.h"
#include "../protocolStack/packet/Packet.h"
#include "../core/eventScheduler/simulator.h"
#include "../flows/application/InfiniteBuffer.h"
#include "../flows/application/VoIP.h"
#include "../flows/application/CBR.h"
#include "../flows/application/TraceBased.h"
#include "../device/IPClassifier/ClassifierParameters.h"
#include "../flows/QoS/QoSParameters.h"
#include "../flows/QoS/QoSForEXP.h"
#include "../flows/QoS/QoSForFLS.h"
#include "../flows/QoS/QoSForM_LWDF.h"
#include "../componentManagers/FrameManager.h"
#include "../utility/seed.h"
#include "../utility/RandomVariable.h"
#include "../phy/wideband-cqi-eesm-error-model.h"
#include "../phy/simple-error-model.h"
#include "../channel/propagation-model/macrocell-urban-area-channel-realization.h"
#include "../load-parameters.h"
#include <queue>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <utility>
using std::pair;
//#define SLICE_CLUSTER 1

static void SingleCellWithInterference (int nbCells, double radius,
                                        int nbUE,
                                        int nbVoIP, int nbVideo, int nbBE, int nbCBR,
                                        int sched_type,
                                        int frame_struct,
                                        int speed,
		                                double maxDelay, int videoBitRate,
                                        int seed,
										string config_fname)
{

  // define simulation times
  double duration = 12;
  double flow_duration = 12;


  int cluster = 3;
  double bandwidth = 20;

  // CREATE COMPONENT MANAGER
  Simulator *simulator = Simulator::Init();
  FrameManager *frameManager = FrameManager::Init();
  NetworkManager* nm = NetworkManager::Init();

  // CONFIGURE SEED
  if (seed >= 0)
	{
	  int commonSeed = GetCommonSeed (seed);
	  srand (commonSeed);
	}
  else
	{
	  srand (time(NULL));
	}
  std::cerr << "Simulation with SEED = " << seed << std::endl;

  // SET SCHEDULING ALLOCATION SCHEME
  ENodeB::DLSchedulerType downlink_scheduler_type;
  switch (sched_type)
	{
	  case 1:
		downlink_scheduler_type = ENodeB::DLScheduler_TYPE_PROPORTIONAL_FAIR;
		std::cout << "Scheduler PF "<< std::endl;
		break;
	  case 2:
		downlink_scheduler_type = ENodeB::DLScheduler_TYPE_MLWDF;
		std::cout << "Scheduler MLWDF "<< std::endl;
		break;
	  case 3:
		downlink_scheduler_type = ENodeB::DLScheduler_TYPE_EXP;
		std::cout << "Scheduler EXP "<< std::endl;
		break;
	  case 4:
		downlink_scheduler_type = ENodeB::DLScheduler_TYPE_FLS;
		std::cout << "Scheduler FLS "<< std::endl;
		break;
    	case 5:
  		downlink_scheduler_type = ENodeB::DLScheduler_EXP_RULE;
  		std::cout << "Scheduler EXP_RULE "<< std::endl;
	  	break;
    	case 6:
  		downlink_scheduler_type = ENodeB::DLScheduler_LOG_RULE;
  		std::cout << "Scheduler LOG RULE "<< std::endl;
	  	break;
    	case 7:
   		downlink_scheduler_type = ENodeB::DLScheduler_NVS;
    	std::cout << "Scheduler NVS " << std::endl;
		break;
		case 8:
		downlink_scheduler_type = ENodeB::DLScheduler_GREEDY;
		std::cout << "Scheduler Greedy " << std::endl;
		break;
		case 9:
		downlink_scheduler_type = ENodeB::DLScheduler_SUBOPT;
		std::cout << "Scheduler Oracle " << std::endl;
		break;
		case 10:
		downlink_scheduler_type = ENodeB::DLScheduler_MAXCELL;
		std::cout << "Scheduler Transport " << std::endl;
		break;
		case 11:
		downlink_scheduler_type = ENodeB::DLScheduler_VOGEL;
		std::cout << "Scheduler VOGEL Approximate";
		break;
	
	  default:
		downlink_scheduler_type = ENodeB::DLScheduler_TYPE_PROPORTIONAL_FAIR;
		break;
	}

  // SET FRAME STRUCTURE
  FrameManager::FrameStructure frame_structure;
  switch (frame_struct)
	{
	  case 1:
		frame_structure = FrameManager::FRAME_STRUCTURE_FDD;
		break;
	  case 2:
		frame_structure = FrameManager::FRAME_STRUCTURE_TDD;
		break;
	  default:
		frame_structure = FrameManager::FRAME_STRUCTURE_FDD;
		break;
	}
  frameManager->SetFrameStructure(FrameManager::FRAME_STRUCTURE_FDD);


  //create cells
  std::vector <Cell*> *cells = new std::vector <Cell*>;
  for (int i = 0; i < nbCells; i++)
	{
	  CartesianCoordinates center =
			  GetCartesianCoordinatesForCell(i, radius * 1000.);

	  Cell *c = new Cell (i, radius, 0.035, center.GetCoordinateX (), center.GetCoordinateY ());
	  cells->push_back (c);
	  nm->GetCellContainer ()->push_back (c);

	  std::cout << "Created Cell, id " << c->GetIdCell ()
			  <<", position: " << c->GetCellCenterPosition ()->GetCoordinateX ()
			  << " " << c->GetCellCenterPosition ()->GetCoordinateY () << std::endl;
	}


  std::vector <BandwidthManager*> spectrums = RunFrequencyReuseTechniques (nbCells, cluster, bandwidth);

  //Create a set of a couple of channels
  std::vector <LteChannel*> *dlChannels = new std::vector <LteChannel*>;
  std::vector <LteChannel*> *ulChannels = new std::vector <LteChannel*>;
  for (int i= 0; i < nbCells; i++)
	{
	  LteChannel *dlCh = new LteChannel ();
	  dlCh->SetChannelId (i);
	  dlChannels->push_back (dlCh);

	  LteChannel *ulCh = new LteChannel ();
	  ulCh->SetChannelId (i);
	  ulChannels->push_back (ulCh);
	}


  //create eNBs
  std::vector <ENodeB*> *eNBs = new std::vector <ENodeB*>;
  for (int i = 0; i < nbCells; i++)
	{
	  ENodeB* enb = new ENodeB (i, cells->at (i));
	  enb->GetPhy ()->SetDlChannel (dlChannels->at (i));
	  enb->GetPhy ()->SetUlChannel (ulChannels->at (i));

	  enb->SetDLScheduler (downlink_scheduler_type, config_fname);

	  enb->GetPhy ()->SetBandwidthManager (spectrums.at (i));

	  std::cout << "Created enb, id " << enb->GetIDNetworkNode()
			  << ", cell id " << enb->GetCell ()->GetIdCell ()
			  <<", position: " << enb->GetMobilityModel ()->GetAbsolutePosition ()->GetCoordinateX ()
			  << " " << enb->GetMobilityModel ()->GetAbsolutePosition ()->GetCoordinateY ()
			  << ", channels id " << enb->GetPhy ()->GetDlChannel ()->GetChannelId ()
			  << enb->GetPhy ()->GetUlChannel ()->GetChannelId ()  << std::endl;

	  spectrums.at (i)->Print ();


	  ulChannels->at (i)->AddDevice((NetworkNode*) enb);


	  nm->GetENodeBContainer ()->push_back (enb);
	  eNBs->push_back (enb);
	}




  //Define Application Container
  int nbCell=1;
  VoIP VoIPApplication[nbVoIP*nbCell*nbUE];
  TraceBased VideoApplication[nbVideo*nbCell*nbUE];
  InfiniteBuffer BEApplication[nbBE*nbCell*nbUE];
  CBR CBRApplication[nbCBR*nbCell*nbUE];
  int voipApplication = 0;
  int videoApplication = 0;
  int cbrApplication = 0;
  int beApplication = 0;
  int destinationPort = 101;
  int applicationID = 0;



  //Create GW
  Gateway *gw = new Gateway ();
  nm->GetGatewayContainer ()->push_back (gw);

  //Create UEs
  int idUE = nbCells;

  #ifdef SLICE_CLUSTER
	int num_slices = 0, ue_id = 0;
	std::vector<int> appid_to_slice;
	std::vector<double> slice_xpos;
	std::vector<double> slice_ypos;
	double temp;
	std::ifstream ifs(config_fname);
	if (!ifs.is_open()) {
		throw std::runtime_error("Defined SLICE_CLUSTER without offering config_file");
	}
	if (ifs.is_open()) {
    	ifs >> num_slices;
    	for (int i = 0; i < num_slices; ++i) {
			ifs >> temp;
	  		double posX = (double)rand()/RAND_MAX;
			double posY = (double)rand()/RAND_MAX;
			int range = 300;
	  		//posX = 0.95 * (((2*radius*1000) * posX) - (radius*1000));
	  		//posY = 0.95 * (((2*radius*1000) * posY) - (radius*1000));
			posX = 0.95 * ((2 * range * radius * posX) - (radius * range) + range * i);
			posY = 0.95 * ((2 * range * radius * posY) - (radius * range) + range * i);
			slice_xpos.push_back(posX);
			slice_ypos.push_back(posY);
		}
		for (int i = 0; i < num_slices; ++i) {
			int ue_num = 0;
			ifs >> ue_num;
			for (int j = 0; j < ue_num; ++j) {
				appid_to_slice.push_back(i);
			}
		}
  }
  ifs.close();
  #endif

  for (int i = 0; i < nbUE; i++)
	{
		#ifdef SLICE_CLUSTER
		int slice_id = appid_to_slice[i];
		//double posX = slice_xpos[slice_id] + (double)rand() / RAND_MAX * radius * 20;
		//double posY = slice_ypos[slice_id] + (double)rand() / RAND_MAX * radius * 20;
		double posX = slice_xpos[slice_id];
		double posY = slice_ypos[slice_id];

		#else
		// double posX = (double)rand()/RAND_MAX;
	  	// posX = 0.95 * (((2*radius*1000)*posX) - (radius*1000));
	  	// double posY = (double)rand()/RAND_MAX;
	  	// posY = 0.95 * (((2*radius*1000)*posY) - (radius*1000));
		double posX = (double)rand() / RAND_MAX * radius * 1000 * 0.6 + 200;
		double posY = (double)rand() / RAND_MAX * radius * 1000 * 0.6 + 200;
		posX = rand() % 2 == 0 ? posX : -posX;
		posY = rand() % 2 == 0 ? posY : -posY;
		#endif
	  //ue's random position
	  

	  double speedDirection = GetRandomVariable (360.) * ((2.*3.14)/360.);

	  UserEquipment* ue = new UserEquipment (idUE,
			                                 posX, posY, speed, speedDirection,
											 cells->at (0),
											 eNBs->at (0),
			                                 0, //handover false!
			                                 //Mobility::CONSTANT_POSITION);
											 Mobility::MANHATTAN);

	  std::cout << "Created UE - id " << idUE << " position " << posX << " " << posY << " direction " << speedDirection << std::endl;

	  ue->GetMobilityModel()->GetAbsolutePosition()->Print();
	  ue->GetPhy ()->SetDlChannel (eNBs->at (0)->GetPhy ()->GetDlChannel ());
	  ue->GetPhy ()->SetUlChannel (eNBs->at (0)->GetPhy ()->GetUlChannel ());

      FullbandCqiManager *cqiManager = new FullbandCqiManager ();
      cqiManager->SetCqiReportingMode (CqiManager::PERIODIC);
      cqiManager->SetReportingInterval (40);
	  //cqiManager->SetReportingInterval (1);
      cqiManager->SetDevice (ue);
      ue->SetCqiManager (cqiManager);

      WidebandCqiEesmErrorModel *errorModel = new WidebandCqiEesmErrorModel ();
      ue->GetPhy ()->SetErrorModel (errorModel);

	  nm->GetUserEquipmentContainer ()->push_back (ue);

	  // register ue to the enb
	  eNBs->at (0)->RegisterUserEquipment (ue);

	  // define the channel realization
	  MacroCellUrbanAreaChannelRealization* c_dl = new MacroCellUrbanAreaChannelRealization (eNBs->at (0), ue);
	  eNBs->at (0)->GetPhy ()->GetDlChannel ()->GetPropagationLossModel ()->AddChannelRealization (c_dl);
	  MacroCellUrbanAreaChannelRealization* c_ul = new MacroCellUrbanAreaChannelRealization (ue, eNBs->at (0));
	  eNBs->at (0)->GetPhy ()->GetUlChannel ()->GetPropagationLossModel ()->AddChannelRealization (c_ul);


	  // CREATE DOWNLINK APPLICATION FOR THIS UE
	  //double start_time = 0.1 + GetRandomVariable (5.);
	  //double duration_time = start_time + flow_duration;
	  double start_time = 0.1;
	  double duration_time = start_time + flow_duration;

	  // *** voip application
	  for (int j = 0; j < nbVoIP; j++)
		{
		  // create application
		  VoIPApplication[voipApplication].SetSource (gw);
		  VoIPApplication[voipApplication].SetDestination (ue);
		  VoIPApplication[voipApplication].SetApplicationID (applicationID);
		  VoIPApplication[voipApplication].SetStartTime(start_time);
		  VoIPApplication[voipApplication].SetStopTime(duration_time);

		  // create qos parameters
		  if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_FLS)
			{
			  QoSForFLS *qos = new QoSForFLS ();
			  qos->SetMaxDelay (maxDelay);
			  if (maxDelay == 0.1)
				{
				  std::cout << "Target Delay = 0.1 s, M = 9" << std::endl;
				  qos->SetNbOfCoefficients (9);
				}
			  else if (maxDelay == 0.08)
				{
				  std::cout << "Target Delay = 0.08 s, M = 7" << std::endl;
				  qos->SetNbOfCoefficients (7);
				}
			  else if (maxDelay == 0.06)
				{
				  std::cout << "Target Delay = 0.06 s, M = 5" << std::endl;
				  qos->SetNbOfCoefficients (5);
				}
			  else if (maxDelay == 0.04)
				{
				  std::cout << "Target Delay = 0.04 s, M = 3" << std::endl;
				  qos->SetNbOfCoefficients (3);
				}
			  else
				{
				  std::cout << "ERROR: target delay is not available"<< std::endl;
				  return;
				}

			  VoIPApplication[voipApplication].SetQoSParameters (qos);
			}
		  else if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_EXP)
			{
			  QoSForEXP *qos = new QoSForEXP ();
			  qos->SetMaxDelay (maxDelay);
			  VoIPApplication[voipApplication].SetQoSParameters (qos);
			}
		  else if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_MLWDF)
			{
			  QoSForM_LWDF *qos = new QoSForM_LWDF ();
			  qos->SetMaxDelay (maxDelay);
			  VoIPApplication[voipApplication].SetQoSParameters (qos);
			}
		  else
			{
			  QoSParameters *qos = new QoSParameters ();
			  qos->SetMaxDelay (maxDelay);
			  VoIPApplication[voipApplication].SetQoSParameters (qos);
			}


		  //create classifier parameters
		  ClassifierParameters *cp = new ClassifierParameters (gw->GetIDNetworkNode(),
															   ue->GetIDNetworkNode(),
															   0,
															   destinationPort,
															   TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
		  VoIPApplication[voipApplication].SetClassifierParameters (cp);

		  std::cout << "CREATED VOIP APPLICATION, ID " << applicationID << std::endl;

		  //update counter
		  destinationPort++;
		  applicationID++;
		  voipApplication++;
		}


	  // *** video application
	  for (int j = 0; j < nbVideo; j++)
		{
		  // create application
		  VideoApplication[videoApplication].SetSource (gw);
		  VideoApplication[videoApplication].SetDestination (ue);
		  VideoApplication[videoApplication].SetApplicationID (applicationID);
		  VideoApplication[videoApplication].SetStartTime(start_time);
		  VideoApplication[videoApplication].SetStopTime(duration_time);

		  string video_trace ("foreman_H264_");
		  //string video_trace ("highway_H264_");
		  //string video_trace ("mobile_H264_");

		  switch (videoBitRate)
			  {
				case 128:
				  {
				    string _file (path + "src/flows/application/Trace/" + video_trace + "128k.dat");
				    VideoApplication[videoApplication].SetTraceFile(_file);
				    std::cout << "		selected video @ 128k " << _file << std::endl;
				    break;
				  }
				case 242:
				  {
				        string _file (path + "src/flows/application/Trace/" + video_trace + "242k.dat");
				    VideoApplication[videoApplication].SetTraceFile(_file);
				    std::cout << "		selected video @ 242k"<< std::endl;
				    break;
				  }
				case 440:
				  {
				    string _file (path + "src/flows/application/Trace/" + video_trace + "440k.dat");
 			        VideoApplication[videoApplication].SetTraceFile(_file);
				    std::cout << "		selected video @ 440k"<< std::endl;
				    break;
				  }
				default:
				  {
				    string _file (path + "src/flows/application/Trace/" + video_trace + "128k.dat");
				    VideoApplication[videoApplication].SetTraceFile(_file);
				    std::cout << "		selected video @ 128k as default"<< std::endl;
				    break;
				  }
			  }

		  // create qos parameters
		  if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_FLS)
			{
			  QoSForFLS *qos = new QoSForFLS ();
			  qos->SetMaxDelay (maxDelay);
			  if (maxDelay == 0.1)
				{
				  std::cout << "Target Delay = 0.1 s, M = 9" << std::endl;
				  qos->SetNbOfCoefficients (9);
				}
			  else if (maxDelay == 0.08)
				{
				  std::cout << "Target Delay = 0.08 s, M = 7" << std::endl;
				  qos->SetNbOfCoefficients (7);
				}
			  else if (maxDelay == 0.06)
				{
				  std::cout << "Target Delay = 0.06 s, M = 5" << std::endl;
				  qos->SetNbOfCoefficients (5);
				}
			  else if (maxDelay == 0.04)
				{
				  std::cout << "Target Delay = 0.04 s, M = 3" << std::endl;
				  qos->SetNbOfCoefficients (3);
				}
			  else
				{
				  std::cout << "ERROR: target delay is not available"<< std::endl;
				  return;
				}
			  VideoApplication[videoApplication].SetQoSParameters (qos);
			}
		  else if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_EXP)
			{
			  QoSForEXP *qos = new QoSForEXP ();
			  qos->SetMaxDelay (maxDelay);
			  VideoApplication[videoApplication].SetQoSParameters (qos);
			}
		  else if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_MLWDF)
			{
			  QoSForM_LWDF *qos = new QoSForM_LWDF ();
			  qos->SetMaxDelay (maxDelay);
			  VideoApplication[videoApplication].SetQoSParameters (qos);
			}
		  else
			{
			  QoSParameters *qos = new QoSParameters ();
			  qos->SetMaxDelay (maxDelay);
			  VideoApplication[videoApplication].SetQoSParameters (qos);
			}


		  //create classifier parameters
		  ClassifierParameters *cp = new ClassifierParameters (gw->GetIDNetworkNode(),
															   ue->GetIDNetworkNode(),
															   0,
															   destinationPort,
															   TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
		  VideoApplication[videoApplication].SetClassifierParameters (cp);

		  std::cout << "CREATED VIDEO APPLICATION, ID " << applicationID << std::endl;

		  //update counter
		  destinationPort++;
		  applicationID++;
		  videoApplication++;
		}

	  // *** be application
	  for (int j = 0; j < nbBE; j++)
		{
		  // create application
		  BEApplication[beApplication].SetSource (gw);
		  BEApplication[beApplication].SetDestination (ue);
		  BEApplication[beApplication].SetApplicationID (applicationID);
		  BEApplication[beApplication].SetStartTime(start_time);
		  BEApplication[beApplication].SetStopTime(duration_time);


		  // create qos parameters
		  QoSParameters *qosParameters = new QoSParameters ();
		  BEApplication[beApplication].SetQoSParameters (qosParameters);


		  //create classifier parameters
		  ClassifierParameters *cp = new ClassifierParameters (gw->GetIDNetworkNode(),
															   ue->GetIDNetworkNode(),
															   0,
															   destinationPort,
															   TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
		  BEApplication[beApplication].SetClassifierParameters (cp);

		  std::cout << "CREATED BE APPLICATION, ID " << applicationID << std::endl;

		  //update counter
		  destinationPort++;
		  applicationID++;
		  beApplication++;
		}

	  // *** cbr application
	  for (int j = 0; j < nbCBR; j++)
		{
		  // create application
		  CBRApplication[cbrApplication].SetSource (gw);
		  CBRApplication[cbrApplication].SetDestination (ue);
		  CBRApplication[cbrApplication].SetApplicationID (applicationID);
		  CBRApplication[cbrApplication].SetStartTime(start_time);
		  CBRApplication[cbrApplication].SetStopTime(duration_time);
		  CBRApplication[cbrApplication].SetInterval (0.04);
		  CBRApplication[cbrApplication].SetSize (5);

		  // create qos parameters
		  QoSParameters *qosParameters = new QoSParameters ();
		  qosParameters->SetMaxDelay (maxDelay);

		  CBRApplication[cbrApplication].SetQoSParameters (qosParameters);


		  //create classifier parameters
		  ClassifierParameters *cp = new ClassifierParameters (gw->GetIDNetworkNode(),
															   ue->GetIDNetworkNode(),
															   0,
															   destinationPort,
															   TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
		  CBRApplication[cbrApplication].SetClassifierParameters (cp);

		  std::cout << "CREATED CBR APPLICATION, ID " << applicationID << std::endl;

		  //update counter
		  destinationPort++;
		  applicationID++;
		  cbrApplication++;
		}

	  idUE++;

	}



  simulator->SetStop(duration);
  simulator->Run ();



  //Delete created objects
  cells->clear ();
  delete cells;
  eNBs->clear ();
  delete eNBs;
  delete frameManager;
  //delete nm;
  delete simulator;

}
