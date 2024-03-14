#!/bin/bash

matlab -nosplash -nodesktop <   localize_tmp/run_1oca_param_1oca_part_1.m &
matlab -nosplash -nodesktop <   localize_tmp/run_1oca_param_1oca_part_2.m &
matlab -nosplash -nodesktop <   localize_tmp/run_1oca_param_1oca_part_3.m &

wait



