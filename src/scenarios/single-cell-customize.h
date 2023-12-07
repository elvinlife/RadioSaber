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

#include <stdlib.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <utility>
#include <vector>

#include "../channel/LteChannel.h"
#include "../channel/propagation-model/macrocell-urban-area-channel-realization.h"
#include "../componentManagers/FrameManager.h"
#include "../core/eventScheduler/simulator.h"
#include "../core/spectrum/bandwidth-manager.h"
#include "../device/IPClassifier/ClassifierParameters.h"
#include "../flows/QoS/QoSForEXP.h"
#include "../flows/QoS/QoSForFLS.h"
#include "../flows/QoS/QoSForM_LWDF.h"
#include "../flows/QoS/QoSParameters.h"
#include "../flows/application/InfiniteBuffer.h"
#include "../flows/application/InternetFlow.h"
#include "../flows/application/TraceBased.h"
#include "../flows/application/VoIP.h"
#include "../load-parameters.h"
#include "../networkTopology/Cell.h"
#include "../phy/enb-lte-phy.h"
#include "../phy/simple-error-model.h"
#include "../phy/ue-lte-phy.h"
#include "../phy/wideband-cqi-eesm-error-model.h"
#include "../protocolStack/packet/Packet.h"
#include "../protocolStack/packet/packet-burst.h"
#include "../utility/RandomVariable.h"
#include "../utility/seed.h"
using std::pair;
using std::vector;

static void SingleCellCustomize(double radius, int sched_type, int frame_struct,
                                int speed, double maxDelay, int videoBitRate,
                                double internetFlowRate, int seed,
                                string config_fname) {
  double duration = 40;
  double flow_duration = 40;

  int nbCells = 1;
  int cluster = 3;
  double bandwidth = 100;

  int nb_be_sliceA = 1;
  int nb_internetflow_sliceB = 1;
  int nb_internetflow_sliceC = 2;
  int nb_videoflow_sliceD = 1;

  // CREATE COMPONENT MANAGER
  Simulator *simulator = Simulator::Init();
  FrameManager *frameManager = FrameManager::Init();
  NetworkManager *nm = NetworkManager::Init();

  // CONFIGURE SEED
  if (seed >= 0) {
    int commonSeed = GetCommonSeed(seed);
    srand(commonSeed);
  } else {
    srand(time(NULL));
  }
  std::cerr << "Simulation with SEED = " << seed << std::endl;

  // SET SCHEDULING ALLOCATION SCHEME
  ENodeB::DLSchedulerType downlink_scheduler_type;
  switch (sched_type) {
    case 1:
      downlink_scheduler_type = ENodeB::DLScheduler_TYPE_PROPORTIONAL_FAIR;
      std::cout << "Scheduler PF " << std::endl;
      break;
    case 2:
      downlink_scheduler_type = ENodeB::DLScheduler_TYPE_MLWDF;
      std::cout << "Scheduler MLWDF " << std::endl;
      break;
    case 3:
      downlink_scheduler_type = ENodeB::DLScheduler_TYPE_EXP;
      std::cout << "Scheduler EXP " << std::endl;
      break;
    case 4:
      downlink_scheduler_type = ENodeB::DLScheduler_TYPE_FLS;
      std::cout << "Scheduler FLS " << std::endl;
      break;
    case 5:
      downlink_scheduler_type = ENodeB::DLScheduler_EXP_RULE;
      std::cout << "Scheduler EXP_RULE " << std::endl;
      break;
    case 6:
      downlink_scheduler_type = ENodeB::DLScheduler_LOG_RULE;
      std::cout << "Scheduler LOG RULE " << std::endl;
      break;
    case 7:
      downlink_scheduler_type = ENodeB::DLScheduler_NVS;
      std::cout << "Scheduler NVS " << std::endl;
      break;
    case 8:
      downlink_scheduler_type = ENodeB::DLScheduler_SEQUENTIAL;
      std::cout << "Scheduler Greedy " << std::endl;
      break;
    case 10:
      downlink_scheduler_type = ENodeB::DLScheduler_UpperBound;
      std::cout << "Scheduler UpperBound " << std::endl;
      break;
    case 11:
      downlink_scheduler_type = ENodeB::DLScheduler_MAXCELL;
      std::cout << "Scheduler MaxCell " << std::endl;
      break;
    case 12:
      downlink_scheduler_type = ENodeB::DLScheduler_VOGEL;
      std::cout << "Scheduler Vogel " << std::endl;
      break;
    default:
      string error_log = "Undefined Scheduler: " + std::to_string(sched_type);
      throw std::runtime_error(error_log);
      break;
  }

  // SET FRAME STRUCTURE
  FrameManager::FrameStructure frame_structure;
  switch (frame_struct) {
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

  // create cells
  std::vector<Cell *> *cells = new std::vector<Cell *>;
  for (int i = 0; i < nbCells; i++) {
    CartesianCoordinates center =
        GetCartesianCoordinatesForCell(i, radius * 1000.);

    Cell *c = new Cell(i, radius, 0.035, center.GetCoordinateX(),
                       center.GetCoordinateY());
    cells->push_back(c);
    nm->GetCellContainer()->push_back(c);

    std::cout << "Created Cell, id " << c->GetIdCell()
              << ", position: " << c->GetCellCenterPosition()->GetCoordinateX()
              << " " << c->GetCellCenterPosition()->GetCoordinateY()
              << std::endl;
  }

  std::vector<BandwidthManager *> spectrums =
      RunFrequencyReuseTechniques(nbCells, cluster, bandwidth);

  // Create a set of a couple of channels
  std::vector<LteChannel *> *dlChannels = new std::vector<LteChannel *>;
  std::vector<LteChannel *> *ulChannels = new std::vector<LteChannel *>;
  for (int i = 0; i < nbCells; i++) {
    LteChannel *dlCh = new LteChannel();
    dlCh->SetChannelId(i);
    dlChannels->push_back(dlCh);

    LteChannel *ulCh = new LteChannel();
    ulCh->SetChannelId(i);
    ulChannels->push_back(ulCh);
  }

  // create eNBs
  std::vector<ENodeB *> *eNBs = new std::vector<ENodeB *>;
  for (int i = 0; i < nbCells; i++) {
    ENodeB *enb = new ENodeB(i, cells->at(i));
    enb->GetPhy()->SetDlChannel(dlChannels->at(i));
    enb->GetPhy()->SetUlChannel(ulChannels->at(i));

    enb->SetDLScheduler(downlink_scheduler_type, config_fname);

    enb->GetPhy()->SetBandwidthManager(spectrums.at(i));

    std::cout
        << "Created enb, id " << enb->GetIDNetworkNode() << ", cell id "
        << enb->GetCell()->GetIdCell() << ", position: "
        << enb->GetMobilityModel()->GetAbsolutePosition()->GetCoordinateX()
        << " "
        << enb->GetMobilityModel()->GetAbsolutePosition()->GetCoordinateY()
        << ", channels id " << enb->GetPhy()->GetDlChannel()->GetChannelId()
        << enb->GetPhy()->GetUlChannel()->GetChannelId() << std::endl;

    spectrums.at(i)->Print();
    ulChannels->at(i)->AddDevice((NetworkNode *)enb);
    nm->GetENodeBContainer()->push_back(enb);
    eNBs->push_back(enb);
  }

  // Define Application Container

  vector<TraceBased *> VideoApplication;
  vector<InfiniteBuffer *> BEApplication;
  vector<InternetFlow *> IPApplication;

  int videoApplication = 0;
  int beApplication = 0;
  int ipApplication = 0;
  int destinationPort = 101;
  int applicationID = 0;

  int num_slices;
  int user_to_slice[1000];
  int slice_users[100];
  int total_ues = 0;
  std::ifstream ifs(config_fname);
  std::string line;
  // read the configuration file
  if (ifs.is_open()) {
    while (std::getline(ifs, line)) {
      if (line[0] == '#') {
        continue;
      }
      std::string args;
      std::istringstream iss(line);
      iss >> args;
      if (args == "num_slices:") {
        iss >> num_slices;
      } else if (args == "slice_ues:") {
        int begin_id = 0;
        for (int i = 0; i < num_slices; ++i) {
          int num_ue;
          iss >> num_ue;
          for (int j = 0; j < num_ue; ++j) {
            user_to_slice[begin_id + j] = i;
          }
          slice_users[i] = num_ue;
          begin_id += num_ue;
          total_ues += num_ue;
        }
      }
    }
  } else {
    throw std::runtime_error("Error, cannot open configuration file.");
  }

  // Create GW
  Gateway *gw = new Gateway();
  nm->GetGatewayContainer()->push_back(gw);

  for (int idUE = 0; idUE < total_ues; idUE++) {
    // we hardcoded the configuration for the customization experiment
    int nbVideo = 0, nbBE = 0, nbInternetFlow = 0;
    bool multi_priority = false;
    int slices_per_group = 5;
    if (user_to_slice[idUE] < slices_per_group) {
      nbBE = nb_be_sliceA;
    } else if (user_to_slice[idUE] < 2 * slices_per_group) {
      nbInternetFlow = nb_internetflow_sliceB;
    } else if (user_to_slice[idUE] < 3 * slices_per_group) {
      nbInternetFlow = nb_internetflow_sliceC;
      multi_priority = true;
    } else if (user_to_slice[idUE] < 4 * slices_per_group) {
      nbVideo = nb_videoflow_sliceD;
    }

    double posX = (double)rand() / RAND_MAX * radius * 1000 * 0.4 + 100;
    double posY = (double)rand() / RAND_MAX * radius * 1000 * 0.4 + 100;
    posX = rand() % 2 == 0 ? posX : -posX;
    posY = rand() % 2 == 0 ? posY : -posY;
    double speedDirection = GetRandomVariable(360.) * ((2. * 3.14) / 360.);

    UserEquipment *ue = new UserEquipment(
        idUE, posX, posY, speed, speedDirection, cells->at(0), eNBs->at(0),
        0,  // handover false!
        Mobility::CONSTANT_POSITION);

    std::cout << "Created UE - id " << idUE << " position " << posX << " "
              << posY << " direction " << speedDirection << std::endl;

    // ue->GetMobilityModel()->GetAbsolutePosition()->Print();
    ue->GetPhy()->SetDlChannel(eNBs->at(0)->GetPhy()->GetDlChannel());
    ue->GetPhy()->SetUlChannel(eNBs->at(0)->GetPhy()->GetUlChannel());

    FullbandCqiManager *cqiManager = new FullbandCqiManager();
    cqiManager->SetCqiReportingMode(CqiManager::PERIODIC);
    cqiManager->SetReportingInterval(40);
    // cqiManager->SetReportingInterval (1);
    cqiManager->SetDevice(ue);
    ue->SetCqiManager(cqiManager);

    WidebandCqiEesmErrorModel *errorModel = new WidebandCqiEesmErrorModel();
    ue->GetPhy()->SetErrorModel(errorModel);

    nm->GetUserEquipmentContainer()->push_back(ue);

    // register ue to the enb
    eNBs->at(0)->RegisterUserEquipment(ue);

    // define the channel realization
    MacroCellUrbanAreaChannelRealization *c_dl =
        new MacroCellUrbanAreaChannelRealization(eNBs->at(0), ue);
    eNBs->at(0)
        ->GetPhy()
        ->GetDlChannel()
        ->GetPropagationLossModel()
        ->AddChannelRealization(c_dl);
    MacroCellUrbanAreaChannelRealization *c_ul =
        new MacroCellUrbanAreaChannelRealization(ue, eNBs->at(0));
    eNBs->at(0)
        ->GetPhy()
        ->GetUlChannel()
        ->GetPropagationLossModel()
        ->AddChannelRealization(c_ul);

    // CREATE DOWNLINK APPLICATION FOR THIS UE
    double start_time = 0.1;
    double duration_time = start_time + flow_duration;

    // *** video application
    for (int j = 0; j < nbVideo; j++) {
      // create application
      TraceBased *video_app = new TraceBased();
      VideoApplication.push_back(video_app);
      video_app->SetSource(gw);
      video_app->SetDestination(ue);
      video_app->SetApplicationID(applicationID);
      video_app->SetStartTime(start_time);
      video_app->SetStopTime(duration_time);

      string video_trace("foreman_H264_");
      // string video_trace ("highway_H264_");
      // string video_trace ("mobile_H264_");

      switch (videoBitRate) {
        case 128: {
          string _file(path + "src/flows/application/Trace/" + video_trace +
                       "128k.dat");
          video_app->SetTraceFile(_file);
          std::cout << "		selected video @ 128k " << _file
                    << std::endl;
          break;
        }
        case 242: {
          string _file(path + "src/flows/application/Trace/" + video_trace +
                       "242k.dat");
          video_app->SetTraceFile(_file);
          std::cout << "		selected video @ 242k" << std::endl;
          break;
        }
        case 440: {
          string _file(path + "src/flows/application/Trace/" + video_trace +
                       "440k.dat");
          video_app->SetTraceFile(_file);
          std::cout << "		selected video @ 440k" << std::endl;
          break;
        }
        case 880: {
          string _file(path + "src/flows/application/Trace/" + video_trace +
                       "880k.dat");
          video_app->SetTraceFile(_file);
          std::cout << "		selected video @ 880k" << std::endl;
          break;
        }
        case 1280: {
          string _file(path + "src/flows/application/Trace/" + video_trace +
                       "1280k.dat");
          video_app->SetTraceFile(_file);
          std::cout << "		selected video @ 1280k" << std::endl;
          break;
        }
        default: {
          string _file(path + "src/flows/application/Trace/" + video_trace +
                       "128k.dat");
          video_app->SetTraceFile(_file);
          std::cout << "		selected video @ 128k as default"
                    << std::endl;
          break;
        }
      }
      // create qos parameters
      if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_FLS) {
        QoSForFLS *qos = new QoSForFLS();
        qos->SetMaxDelay(maxDelay);
        if (maxDelay == 0.1) {
          std::cout << "Target Delay = 0.1 s, M = 9" << std::endl;
          qos->SetNbOfCoefficients(9);
        } else if (maxDelay == 0.08) {
          std::cout << "Target Delay = 0.08 s, M = 7" << std::endl;
          qos->SetNbOfCoefficients(7);
        } else if (maxDelay == 0.06) {
          std::cout << "Target Delay = 0.06 s, M = 5" << std::endl;
          qos->SetNbOfCoefficients(5);
        } else if (maxDelay == 0.04) {
          std::cout << "Target Delay = 0.04 s, M = 3" << std::endl;
          qos->SetNbOfCoefficients(3);
        } else {
          std::cout << "ERROR: target delay is not available" << std::endl;
          return;
        }
        video_app->SetQoSParameters(qos);
      } else if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_EXP) {
        QoSForEXP *qos = new QoSForEXP();
        qos->SetMaxDelay(maxDelay);
        video_app->SetQoSParameters(qos);
      } else if (downlink_scheduler_type == ENodeB::DLScheduler_TYPE_MLWDF) {
        QoSForM_LWDF *qos = new QoSForM_LWDF();
        qos->SetMaxDelay(maxDelay);
        video_app->SetQoSParameters(qos);
      } else {
        QoSParameters *qos = new QoSParameters();
        qos->SetMaxDelay(maxDelay);
        video_app->SetQoSParameters(qos);
      }

      // create classifier parameters
      ClassifierParameters *cp = new ClassifierParameters(
          gw->GetIDNetworkNode(), ue->GetIDNetworkNode(), 0, destinationPort,
          TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
      video_app->SetClassifierParameters(cp);

      std::cout << "CREATED VIDEO APPLICATION, ID " << applicationID
                << std::endl;

      // update counter
      destinationPort++;
      applicationID++;
      videoApplication++;
    }

    // *** backlogged application
    for (int j = 0; j < nbBE; j++) {
      // create application
      InfiniteBuffer *be_app = new InfiniteBuffer();
      BEApplication.push_back(be_app);
      be_app->SetSource(gw);
      be_app->SetDestination(ue);
      be_app->SetApplicationID(applicationID);
      be_app->SetStartTime(start_time);
      be_app->SetStopTime(duration_time);

      // create qos parameters
      QoSParameters *qosParameters = new QoSParameters();
      be_app->SetQoSParameters(qosParameters);

      // create classifier parameters
      ClassifierParameters *cp = new ClassifierParameters(
          gw->GetIDNetworkNode(), ue->GetIDNetworkNode(), 0, destinationPort,
          TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
      be_app->SetClassifierParameters(cp);

      std::cout << "CREATED BE APPLICATION, ID " << applicationID << std::endl;

      // update counter
      destinationPort++;
      applicationID++;
      beApplication++;
    }

    // *** constant bitrate flows following heavy-tail distributions
    for (int j = 0; j < nbInternetFlow; j++) {
      double user_rate = internetFlowRate / (slice_users[user_to_slice[idUE]]);
      InternetFlow *ip_app = new InternetFlow();
      IPApplication.push_back(ip_app);
      ip_app->SetSource(gw);
      ip_app->SetDestination(ue);
      ip_app->SetApplicationID(applicationID);
      ip_app->SetStartTime(start_time);
      ip_app->SetStopTime(duration_time);
      if (multi_priority) {
        ip_app->SetPriority(j);
        if (j == 0) {
          ip_app->SetAvgRate(user_rate * 0.75);
        } else if (j == 1) {
          ip_app->SetAvgRate(user_rate * 0.25);
        }
      }
      // evenly divide average flow rate
      else {
        ip_app->SetPriority(0);
        ip_app->SetAvgRate(user_rate / nbInternetFlow);
      }

      QoSParameters *qosParameters = new QoSParameters();
      ip_app->SetQoSParameters(qosParameters);

      ClassifierParameters *cp = new ClassifierParameters(
          gw->GetIDNetworkNode(), ue->GetIDNetworkNode(), 0, destinationPort,
          TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
      ip_app->SetClassifierParameters(cp);

      std::cout << "CREATED InternetFlow APPLICATION, ID " << applicationID
                << std::endl;

      destinationPort++;
      applicationID++;
      ipApplication++;
    }
  }

  simulator->SetStop(duration);
  simulator->Run();

  // Delete created objects
  cells->clear();
  delete cells;
  eNBs->clear();
  delete eNBs;
  delete frameManager;
  // delete nm;
  delete simulator;
  delete dlChannels;
  delete ulChannels;

  for (auto it = VideoApplication.begin(); it != VideoApplication.end(); ++it)
    delete *it;
  for (auto it = BEApplication.begin(); it != BEApplication.end(); ++it)
    delete *it;
  for (auto it = IPApplication.begin(); it != IPApplication.end(); ++it)
    delete *it;
}
