RadioSaber is State-of-the-Art 5G RAN slicing algorithm that achieves high spectrum efficiency, ensures weighted fairness among slices subject to their SLA(service-level-agreement), and allow slices to customize their scheduling policies. For more details, please check our paper published in NSDI'2023.

## Overview of RadioSaber
RadioSaber is built upon LTE-Sim, an open source framework to simulate the Radio Access networks in LTE. LTE-Sim simulates all the software stacks in RAN(PDCP, RLC, MAC and RRC), channel propagation model in the physical layer, diverse applications with different QoS requirements, user mobility, etc. Here, RadioSaber mostly focuses on the flow scheduling policy in the MAC layer. To capture the flow scheduling scenario of 5G base stations in real world, we extend the LTE-Sim to support 100MHz downlink in cells, allocation of RBGs(resource block group) instead of RBs(resource block), subband CQI report mechanism, and to use real channel quality traces collected with SDR(software defined radios)

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
  * 11: NVS-Nongreedy, the inter-slice scheduler applies NVS while enterprise schedulers apply non-greedy proportional fairness algorithm proposed by [Mobicom18](https://dl.acm.org/doi/abs/10.1145/3241539.3241552)
* frame struct: the cellular applies FDD or TDD. It must be set to 1 to apply FDD
* mobility speed: the mobility speed of UEs. It doesn't matter here since the CQI of UEs are traced collected from SDR instead of being simulated by the simulator
* random seed: random seed
* config file: a json-based configuration file

Config File:

The config file is in JSON file format. It configures the scheduling algorithm, number of UEs, weight of every slice, and how UEs in every slice instantiate applications and flows.

## Experiments and Reproducibility

### Spectrum Efficiency and Fairness($\S6.1$)
### Diverse Enterprise Schedulers($\S6.2$)
### What Makes RadioSaber Win Over NVS($\S6.3$)
### Varying Number of Slices and UEs per Slice($\S6.4$)
### Non-greedy Enterprise Schedulers($\S6.5$)
### Is There a Better Inter-Slice Scheduler?($\S6.6$)
### Scheduling Latency and Overhead($\S6.7$)

[https://github.com/elvinlife/radiosaber-overhead](https://github.com/elvinlife/radiosaber-overhead)
