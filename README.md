<!--
 * @Author       : Jie Wu j.wu@cern.ch
 * @Date         : 2024-07-19 14:13:59 +0200
 * @LastEditors  : Jie Wu j.wu@cern.ch
 * @LastEditTime : 2024-11-06 08:26:47 +0100
 * @FilePath     : README.md
 * @Description  : 
 * 
 * Copyright (c) 2024 by everyone, All Rights Reserved. 
-->
* PIDCalibration wrapper.

Use correct samples to do PID calibration automatically according to the TRUEID of the particles, based on the tool described in `https://twiki.cern.ch/twiki/bin/viewauth/LHCb/MeerkatPIDResampling`

The following branches should exists in the tree:
1. [track name]_[TRUEID]
2. [nTracks]
3. [track name]_[P]
4. [track name]_[PT]



folders:
- `config` : configuration files for PIDCalib
- `scripts` : core scripts for PIDCalib



scripts:
- `scripts/PIDCalib.py` : Main script for PIDCalib
- `scripts/PIDCorr_git_edit.py` : The core script for PIDCalib (do not need to be modified)


config:
- `config/config_samples_PIDCorr.yaml` : configuration file for PID correction (transfomration)
- `config/config_samples_PIDGen.yaml` : configuration file for PID generation (resampling)
- `config/pidcorr_run1.yaml` : The file which defines the variables you want to correct for run2 MC
- `config/pidcorr_run2.yaml` : The file which defines the variables you want to correct for run2 MC



how to run?
- Source to enter lhcb environment by calling: source /cvmfs/lhcb.cern.ch/lib/LbEnv

- PIDCorr:  sh run_pidcorr_local.sh
- PIDGen:   sh run_pidgen_local.sh



what should I modify?
- The files in the config folder
- Modify `PIDCorr` and `PIDGen` to your needs




The complete list of PID configurations available in Run 1 and Run 2 can be obtained by running PIDGen.py without arguments:
$ lb-run -c best --allow-containers --siteroot=/cvmfs/lhcb.cern.ch/lib Urania/v10r1 bash
$ python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDGen.py
Note that Run 1 and Run 2 PID configurations have different names! 