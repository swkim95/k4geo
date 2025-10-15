#!/bin/bash

cd /u/user/swkim/KEY4HEP/2025/for_intern_recipe/k4geo
source setenv.sh
cd example

ddsim --steeringFile SteeringFile_INCREADIBLE_v01_20GeV_elec.py --outputFile ./test_INCREADIBLE_o1_v01_500Evt_20GeV_electron_${1}.root
