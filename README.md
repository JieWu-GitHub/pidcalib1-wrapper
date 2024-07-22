* PIDCalibration wrapper.

Use correct samples to do PID calibration automatically according to the TRUEID of the particles.

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
- PIDCorr:  sh run_pidcorr_local.sh
- PIDGen:   sh run_pidgen_local.sh



what should I modify?
- The files in the config folder
- Modify `PIDCorr` and `PIDGen` to your needs


