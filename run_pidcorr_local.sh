#!/bin/bash
###
 # @Author       : Jie Wu j.wu@cern.ch
 # @Date         : 2024-07-22 05:34:04 +0200
 # @LastEditors  : Jie Wu j.wu@cern.ch
 # @LastEditTime : 2024-07-23 14:46:53 +0200
 # @FilePath     : run_pidcorr_local.sh
 # @Description  : 
 # 
 # Copyright (c) 2024 by everyone, All Rights Reserved. 
### 




mode=MC_Bs2JpsiKstar
year=2018
pol=MagDown



input_file=test/test_Bd2JpsiKstar_${year}_${pol}.root
input_tree_name=DecayTree



# resamp_var=Pi_PIDmu
resamp_var=${mode}
output_file="test/output/PIDCorr/${year}/${pol}/${mode}.root" 

output_var_suffix="PIDCorr"


# To determine the yaml with resample names 
if [[ $year == "2015" || $year == "2016" || $year == "2017" || $year == "2018" ]]; then
    run="run2"
    calibconfig="run2"
elif [[ $year == "2011" || $year == "2012" ]]; then
    run="run1"
    calibconfig="sim09"
fi

tracks_file="config/PIDCorr/pidcorr_${run}.yaml"
config_file="config/PIDCorr/config_samples_PIDCorr.yaml"


echo "Using yaml file for run: ${run}"



################# COMMAND #################
lb-run --bind=/DATA:/DATA --bind=/EOS:/EOS -c best --platform=x86_64_v2-centos7-gcc11-opt --siteroot=/cvmfs/lhcb.cern.ch/lib Urania/v10r1 python -u scripts/PIDCorr.py \
--input-file "${input_file}" \
--input-tree-name "${input_tree_name}"  \
--output-file  "${output_file}"   \
--data-set ${pol}_${year}           \
--simversion ${calibconfig}         \
--mode ${resamp_var}                  \
--tracks-file ${tracks_file}        \
--config-file ${config_file}        \
--tmp1 "test/output/PIDCorr/tmpc/${year}/${pol}/${mode}_tmp1.root"  \
--tmp2 "test/output/PIDCorr/tmpc/${year}/${pol}/${mode}_tmp2.root"  \
--output-var-suffix ${output_var_suffix}








