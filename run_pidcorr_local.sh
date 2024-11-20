#!/bin/bash
###
 # @Author       : Jie Wu j.wu@cern.ch
 # @Date         : 2024-07-22 05:34:04 +0200
 # @LastEditors  : Jie Wu j.wu@cern.ch
 # @LastEditTime : 2024-11-07 07:41:43 +0100
 # @FilePath     : run_pidcorr_local.sh
 # @Description  : 
 # 
 # Copyright (c) 2024 by everyone, All Rights Reserved. 
### 




mode=MC_Bs2JpsiKstar
year=2018
pol=MagUp

charge_suffix=""
# charge_suffix="_P"
# charge_suffix="_M"


input_file=test/test_Bd2JpsiKstar_${year}_${pol}.root
# input_file=test/match_ppimumu_${year}_${pol}_sizeReduced.root
input_tree_name=DecayTree


# resamp_var=Pi_PIDmu
resamp_var=${mode}
output_file="test/output/PIDCorr/${year}/${pol}/${mode}${charge_suffix}.root" 

output_var_suffix="PIDCorr"


# To determine the yaml with resample names 
if [[ $year == "2015" || $year == "2016" || $year == "2017" || $year == "2018" ]]; then
    calibconfig="run2"
elif [[ $year == "2011" || $year == "2012" ]]; then
    calibconfig="sim09" # sim08, sim09 for run1
fi

tracks_file="config/PIDCorr_PlusMinus/pidcorr_${calibconfig}.yaml"
config_file="config/PIDCorr_PlusMinus/config_samples_PIDCorr.yaml"

local_root_dir="/disk/lhcb_data/jwu/Bs2JpsiKstar/PIDCalibSamples/eos/lhcb/wg/PID/PIDGen"
# local_root_dir="/home/uzh/wjie/workspace/Bs2JpsiKst-fullRun2/PID_asymmetry/PIDCalib/PIDPerfScripts/python/output/eos/lhcb/wg/PID/PIDGen"
# local_root_dir="/EOS/lhcb/wg/PID/PIDGen"
# local_root_dir="NONE"

echo "Using yaml file for calibconfig: ${calibconfig}"



################# COMMAND #################
# lb-run --bind=/DATA:/DATA --bind=/EOS:/EOS -c best --platform=x86_64_v2-centos7-gcc11-opt --siteroot=/cvmfs/lhcb.cern.ch/lib Urania/v10r1 python -u scripts/PIDCorr.py \
# lb-run --bind=/disk:/disk -c best --allow-containers --siteroot=/cvmfs/lhcb.cern.ch/lib Urania/v10r1 python -u scripts/PIDCorr.py \
# lb-run --bind=/disk/lhcb_data/jwu:/disk/lhcb_data/jwu -c best --allow-containers --siteroot=/cvmfs/lhcb.cern.ch/lib Urania/v10r1 python -u scripts/PIDCorr.py \


lb-run --bind=/disk/lhcb_data/jwu:/disk/lhcb_data/jwu -c best --allow-containers --siteroot=/cvmfs/lhcb.cern.ch/lib Urania/v10r1 python -u scripts/PIDCorr.py \
--input-file        "${input_file}" \
--input-tree-name   "${input_tree_name}"  \
--output-file       "${output_file}"   \
--output-tree       "DecayTree"  \
--data-set      "${pol}_${year}${charge_suffix}"           \
--simversion    "${calibconfig}"         \
--mode          "${resamp_var}"                  \
--tracks-file   "${tracks_file}"        \
--config-file   "${config_file}"        \
--tmp1  "test/output/PIDCorr/tmpc/${year}/${pol}/${mode}${charge_suffix}_tmp1.root"  \
--tmp2  "test/output/PIDCorr/tmpc/${year}/${pol}/${mode}${charge_suffix}_tmp2.root"  \
--output-var-suffix "${output_var_suffix}" \
--local-root-dir    "${local_root_dir}"








