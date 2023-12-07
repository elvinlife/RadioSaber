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

#include <jsoncpp/json/json.h>
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

struct SliceConfig {
  int nb_video;
  int nb_internetflow;
  int nb_backlogflow;
  std::vector<int> video_bitrate;  // the video bitrate of a single user
  std::vector<double> if_bitrate;  // the aggregate data rate of the slice

  SliceConfig(int _nb_video, int _nb_internetflow, int _nb_backlogflow)
      : nb_video(_nb_video),
        nb_internetflow(_nb_internetflow),
        nb_backlogflow(_nb_backlogflow) {}
};

static void SingleCellWithInterference(double radius, int sched_type,
                                       int frame_struct, int speed, int seed,
                                       double duration, string config_fname) {
  int nbCells = 1;
  int cluster = 3;
  double bandwidth = 100;

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
    case 7:
      downlink_scheduler_type = ENodeB::DLScheduler_NVS;
      std::cout << "Scheduler NVS " << std::endl;
      break;
    case 8:
      downlink_scheduler_type = ENodeB::DLScheduler_SEQUENTIAL;
      std::cout << "Scheduler Greedy " << std::endl;
      break;
    case 9:
      downlink_scheduler_type = ENodeB::DLScheduler_MAXCELL;
      std::cout << "Scheduler MaxCell " << std::endl;
      break;
    case 10:
      downlink_scheduler_type = ENodeB::DLScheduler_UpperBound;
      std::cout << "Scheduler UpperBound " << std::endl;
      break;
    case 11:
      downlink_scheduler_type = ENodeB::DLScheduler_NVS_NONGREEDY;
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

  int destinationPort = 101;
  int applicationID = 0;
  int num_slices = 0, num_ue = 0;

  vector<int> user_to_slice;
  vector<int> slice_users;
  vector<SliceConfig> slice_configs;
  int total_ues = 0;

  std::ifstream ifs(config_fname);
  Json::Reader reader;
  Json::Value obj;
  if (!reader.parse(ifs, obj)) {
    throw std::runtime_error("Error, failed to parse the json config file.");
  }
  const Json::Value &ues_per_slice = obj["ues_per_slice"];
  num_slices = ues_per_slice.size();
  std::cout << "num_slices:" << num_slices << std::endl;
  for (int i = 0; i < num_slices; i++) {
    num_ue = ues_per_slice[i].asInt();
    for (int j = 0; j < num_ue; j++) {
      user_to_slice.push_back(i);
    }
    slice_users.push_back(num_ue);
    total_ues += num_ue;
  }
  const Json::Value &slice_schemes = obj["slices"];
  for (int i = 0; i < slice_schemes.size(); i++) {
    int n_slices = slice_schemes[i]["n_slices"].asInt();
    for (int j = 0; j < n_slices; j++) {
      SliceConfig config(slice_schemes[i]["video_app"].asInt(),
                         slice_schemes[i]["internet_flow"].asInt(),
                         slice_schemes[i]["backlog_flow"].asInt());
      const Json::Value &video_bitrate = slice_schemes[i]["video_bitrate"];
      for (int k = 0; k < video_bitrate.size(); k++) {
        config.video_bitrate.push_back(video_bitrate[k].asInt());
      }
      const Json::Value &if_bitrate = slice_schemes[i]["if_bitrate"];
      for (int k = 0; k < if_bitrate.size(); k++) {
        config.if_bitrate.push_back(if_bitrate[k].asDouble());
      }
      slice_configs.push_back(config);
    }
  }

  // Create GW
  Gateway *gw = new Gateway();
  nm->GetGatewayContainer()->push_back(gw);

  double start_time = 0.1;
  double duration_time = start_time + duration;

  for (int idUE = 0; idUE < total_ues; idUE++) {
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

    ue->GetMobilityModel()->GetAbsolutePosition()->Print();
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
    SliceConfig config = slice_configs[user_to_slice[idUE]];
    // *** video application
    for (int j = 0; j < config.nb_video; j++) {
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

      switch (config.video_bitrate[j]) {
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
      // create classifier parameters
      ClassifierParameters *cp = new ClassifierParameters(
          gw->GetIDNetworkNode(), ue->GetIDNetworkNode(), 0, destinationPort,
          TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
      video_app->SetClassifierParameters(cp);

      std::cout << "CREATED VIDEO APPLICATION, ID " << applicationID
                << " Video Bitrate: " << config.video_bitrate[j] << " Kbps"
                << std::endl;

      // update counter
      destinationPort++;
      applicationID++;
    }

    // *** backlogged application
    for (int j = 0; j < config.nb_backlogflow; j++) {
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
    }

    // heavy tail distributed internet flows with const bitrate
    for (int j = 0; j < config.nb_internetflow; j++) {
      double flow_rate =
          config.if_bitrate[j] / slice_users[user_to_slice[idUE]];
      InternetFlow *ip_app = new InternetFlow();
      IPApplication.push_back(ip_app);
      ip_app->SetSource(gw);
      ip_app->SetDestination(ue);
      ip_app->SetApplicationID(applicationID);
      ip_app->SetStartTime(start_time);
      ip_app->SetStopTime(duration_time);
      ip_app->SetAvgRate(flow_rate);
      ip_app->SetPriority(j);

      QoSParameters *qosParameters = new QoSParameters();
      ip_app->SetQoSParameters(qosParameters);

      ClassifierParameters *cp = new ClassifierParameters(
          gw->GetIDNetworkNode(), ue->GetIDNetworkNode(), 0, destinationPort,
          TransportProtocol::TRANSPORT_PROTOCOL_TYPE_UDP);
      ip_app->SetClassifierParameters(cp);

      std::cout << "CREATED InternetFlow APPLICATION, ID " << applicationID
                << " Flow Rate: " << flow_rate << " Mbps" << std::endl;

      destinationPort++;
      applicationID++;
    }
  }

  simulator->SetStop(duration_time);
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
