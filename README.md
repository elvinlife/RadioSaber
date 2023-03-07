# RadioSaber
![Status](https://img.shields.io/badge/Version-Experimental-green.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[Website](https://radiosaber.web.illinois.edu/) [Paper](https://www.usenix.org/conference/nsdi23/presentation/chen-yongzhou)

RadioSaber is State-of-the-Art 5G RAN slicing algorithm that achieves high spectrum efficiency, ensures weighted fairness among slices subject to their SLA(service-level-agreement), and allows slices to customize their scheduling policies. For more details, please check our paper published in NSDI'2023.

- [RadioSaber](#radiosaber)
  * [Overview of RadioSaber](#overview-of-radiosaber)
  * [Software Installation](#software-installation)
  * [How to run RadioSaber](#how-to-run-radiosaber)
  * [Experiments and Reproducibility](#experiments-and-reproducibility)
    + [Spectrum Efficiency and Fairness(Sec 6.1):](#spectrum-efficiency-and-fairness-sec-61--)
    + [Diverse Enterprise Schedulers(Sec 6.2):](#diverse-enterprise-schedulers-sec-62--)
    + [What Makes RadioSaber Win Over NVS(Sec 6.3):](#what-makes-radiosaber-win-over-nvs-sec-63--)
    + [Varying Number of Slices and UEs per Slice(Sec 6.4)](#varying-number-of-slices-and-ues-per-slice-sec-64-)
    + [Non-greedy Enterprise Schedulers(Sec 6.5)](#non-greedy-enterprise-schedulers-sec-65-)
    + [Is There a Better Inter-Slice Scheduler?(Sec 6.6)](#is-there-a-better-inter-slice-scheduler--sec-66-)
    + [Scheduling Latency and Overhead(Sec 6.7)](#scheduling-latency-and-overhead-sec-67-)
  * [Contact](#contact)

## Overview of RadioSaber
RadioSaber is built upon LTE-Sim, an open-source framework to simulate the Radio Access networks in LTE. LTE-Sim simulates all the software stacks in RAN(PDCP, RLC, MAC, and RRC), channel propagation model in the physical layer, diverse applications with different QoS requirements, user mobility, etc. Here, RadioSaber mostly focuses on the flow scheduling policy in the MAC layer. To capture the flow scheduling scenario of 5G base stations in the real world, we extend the LTE-Sim to support 100MHz downlink in cells, allocation of RBGs(resource block group) instead of RBs(resource block), subband CQI report mechanism, and use real channel quality traces collected with SDR(software defined radios)

## Software Installation
Currently, RadioSaber is only tested and supported in Ubuntu environment

```
sudo apt-get install libjsoncpp-dev
git clone https://github.com/elvinlife/RadioSaber
cd RadioSaber
make -j8
```

## How to run RadioSaber
Since RadioSaber is built upon LTE-Sim, there is a lot legacy code we won't use. After you've built RadioSaber, run the following command to start an experiment:

```
${PATH-TO-RADIOSABER}/LTE-Sim SingleCellWithI [radius] [inter-slice scheduler] [frame struct] [mobility speed] [random seed] [config file]
```

Parameters:

* radius: the radius of the coverage area of a cell, the default value is 1
* inter-slice scheduler:
  * 1: No-Slicing, a global proportional fairness scheduler without any notion of slicing
  * 7: NVS, a channel-agnostic inter-slice scheduler proposed in [Mobicom10](https://dl.acm.org/doi/10.1145/1859995.1860023)
  * 8: Sequential, a channel-aware inter-slice scheduler that has lower time complexity and lower spectrum efficiency
  * 9: RadioSaber, our channel-aware inter-slice scheduler
  * 10: Upperbound, an impractical scheme that offers an upper-bound on the spectrum efficiency that any inter-slice scheduler can achieve
  * 11: NVS-Nongreedy, the inter-slice scheduler applies NVS while enterprise schedulers apply a non-greedy proportional fairness algorithm proposed by [Mobicom18](https://dl.acm.org/doi/abs/10.1145/3241539.3241552)
* frame struct: the cellular applies FDD or TDD. It must be set to 1 to apply FDD
* mobility speed: the mobility speed of UEs. If the simulator uses collected CQI traces instead of simulating the channel quality, then it doesn't matter.
* random seed: random seed
* config file: a JSON-based configuration file

Config File:

The config file is in JSON file format. It configures the scheduling algorithm, number of UEs, the weight of every slice, and how UEs in every slice instantiate applications and flows.

## Experiments and Reproducibility
### Spectrum Efficiency and Fairness(Sec 6.1):
* Change the working directory: ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/exp-customization```, and run this experiment with ```./run_backlogged.sh```

* For slices with the same weights, the config file is: ```./exp-backlogged-20slicesdiffw/config.json```
For slices with different weights, the config file is: ```./exp-backlogged-20slicesdiffw/config ```

* After experiments finish, run ```./plot_throughput.py``` to get the throughput and radio resource distribution graph

### Diverse Enterprise Schedulers(Sec 6.2):
* Change the working directory: ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/exp-customization```, and run this experiment with ```./run_customize.sh```
* The config file locates in ```./exp-customize-20slices/config.json```
* Run ```./plot_fctdelay.py``` to get the CDF graphs of flow completion time and packet delay time

### What Makes RadioSaber Win Over NVS(Sec 6.3):
For both synthetic experiments, we don't use real CQI traces but leverage the simulator to simulate the channel quality of every user based on his mobility and geo-location.

* Check the global config file ```./src/load-parameters.h```, and comment the flag ```#define USE_REAL_TRACE```. To run the first experiment, uncomment and define the flag ```#define FIRST_SYNTHETIC_EXP```; To run the second experiment, uncomment and define the flag ```#define SECOND_SYNTHETIC_EXP```
* Change the working directory: ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/wideband-100500-urban```, and run ```./run_backlogged.sh``` to get results of the first synthetic experiment
* Change the working directory: ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/noncomplement-100500-urban```, and run ```./run_backlogged.sh``` to get results of the second synthetic experiment

### Varying Number of Slices and UEs per Slice(Sec 6.4)
* Vary the number of slices with 5-15 users per slice:
	* Change the working directory: ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/exp-fixranues```
	* Run the script: ```./run_exps.sh```
* Vary the number of users per slice with 20 slices:
	* Change the working directory: ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/exp-fix20slices```
	* Run the script: ```./run_exps.sh```

### Non-greedy Enterprise Schedulers(Sec 6.5)
The scripts and config files locate in the directory ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/exp-nongreedy```

### Is There a Better Inter-Slice Scheduler?(Sec 6.6)
* Change the working directory: ```cd ${PATH-TO-RADIOSABER}/NSDI23-radiosaber-experiments/exp-fixranues```
* Launch the experiment with: ```./run_multi_interslice.sh```
* Use the function ```plot_exp2_graph()``` in ```plot_all.py``` to plot the graph

### Scheduling Latency and Overhead(Sec 6.7)
[https://github.com/elvinlife/radiosaber-overhead](https://github.com/elvinlife/radiosaber-overhead)

## Contact
Contact yc28 [at] illinois [dot] edu for assistance
