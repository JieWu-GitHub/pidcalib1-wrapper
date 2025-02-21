'''
Author       : Jie Wu j.wu@cern.ch
Date         : 2024-07-22 05:33:10 +0200
LastEditors  : Jie Wu j.wu@cern.ch
LastEditTime : 2024-11-21 03:07:37 +0100
FilePath     : PIDCorr.py
Description  : 

Copyright (c) 2024 by everyone, All Rights Reserved. 
'''

import os
import sys
import yaml
import argparse
from pprint import pprint


def read_from_yaml(mode, selection_files):
    selection_dict = {}
    for f in selection_files:
        print(f"Reading from {f}")
        with open(f, 'r') as stream:
            selection_dict.update(yaml.safe_load(stream)[mode])
    return selection_dict


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', help='Path to the input file')
    parser.add_argument('--input-tree-name', help='Name of the tree')
    parser.add_argument('--output-file', help='Output ROOT file')
    parser.add_argument('--output-tree', help='Name of the output tree')
    parser.add_argument('--data-set', help='Mag and Year, for example MagUp_2016')
    parser.add_argument('--mode', help='Name of the selection in yaml')
    parser.add_argument('--simversion', help='Simversion, where to read samples from')
    parser.add_argument('--tracks-file', nargs='+', help='Yaml file with tracks')
    parser.add_argument('--config-file', help='Path to the PID config sample definition file')
    parser.add_argument('--tmp1', help='Temporary file to apply PID 1st step')
    parser.add_argument('--tmp2', help='Temporary file to apply PID 2nd step')
    parser.add_argument(
        '--output-var-suffix',
        default='PIDCalib',
        help='Suffix for output PIDCalib variables, like if PIDCalib is the suffix, the output variable will be like "mu_PIDmu_PIDCalib" for muon PIDmu calibration',
    )
    parser.add_argument(
        '--local-root-dir',
        default=None,
        help='Path to the local directory where the ROOT files are stored',
    )

    return parser


def PIDCorr(
    input_file,
    input_tree_name,
    output_file,
    output_tree,
    data_set,
    mode,
    simversion,
    tracks_file,
    config_file,
    tmp1,
    tmp2,
    output_var_suffix,
    local_root_dir,
):
    ## START OF CONFIG
    # Read comments and check vars
    # at least until end of config section

    # List of input ROOT files with MC ntuples. Format:
    #   (inputfile, outputfile, dataset)
    files = [
        (input_file, output_file, data_set),
    ]

    # Name of the input tree
    # Could also include ROOT directory, e.g. "Dir/Ntuple"
    input_tree = input_tree_name

    # Set to config from run_PIDCorrection script CALIBCONFIG variable
    simversion = simversion

    # Postfixes of the Pt, Eta and Ntracks variables (ntuple variable name w/o branch name)
    # e.g. if the ntuple contains "pion_PT", it should be just "PT"
    ptvar = "PT"
    etavar = None
    pvar = "P"
    ## Could use P variable instead of eta
    # etavar = None
    # pvar   = "p"

    ntrvar = "nTracks"  # This should correspond to the number of "Best tracks", not "Long tracks"!

    # Dictionary of tracks with their PID variables, in the form {branch name}:{pidvars}
    # For each track branch name, {pidvars} is a dictionary in the form {ntuple variable}:{pid config},
    #   where
    #     {ntuple variable} is the name of the corresponding ntuple PID variable without branch name,
    #   and
    #     {pid_config} is the string describing the PID configuration.
    # Run PIDCorr.py without arguments to get the full list of PID configs
    tracks = read_from_yaml(mode, tracks_file)
    print("Tracks:")
    pprint(tracks)

    # IF ON LXPLUS: if /tmp exists and is accessible, use for faster processing
    # IF NOT: use /tmp if you have enough RAM
    # temp_folder = '/tmp'
    # ELSE: use current folder

    ## END OF CONFIG

    # make sure we don't overwrite local files and prefix them with random strings
    import string
    import random

    rand_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))  # get 10 random chars for temp_file prefix

    _input_tree = input_tree

    for input_file, output_file, dataset in files:
        tmpinfile = input_file
        tmpoutfile = tmp1
        for track, subst in tracks.items():
            # for var, config in subst.items():
            for var, oldpidvar in subst.items():
                oldpidvar = var if oldpidvar == '' else oldpidvar

                # command = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDCorr.py"
                # command = "python PIDCorr_git.py" #Copied from git, to make local edits
                command = "python scripts/PIDCorr_edit.py"  # Copied from git, to make local edits
                command += " -m %s_%s" % (track, ptvar)
                if etavar:
                    command += " -e %s_%s" % (track, etavar)
                elif pvar:
                    command += " -q %s_%s" % (track, pvar)
                else:
                    print('Specify either ETA or P branch name per track')
                    sys.exit(1)
                command += " --ntrvar %s" % ntrvar
                command += " --tree %s" % _input_tree
                command += " --pidvar %s_%s_%s" % (track, var, output_var_suffix)
                command += " --simpidvar %s_%s" % (track, oldpidvar)
                command += " --config %s" % config_file
                command += " --dataset %s" % dataset
                command += " --input %s" % tmpinfile
                command += " --output %s" % tmpoutfile
                command += " --simversion %s" % simversion
                command += " --outtree %s" % output_tree
                command += " --localrootdir %s" % local_root_dir

                _input_tree = output_tree
                tmpinfile = tmpoutfile
                if tmpoutfile == tmp1:
                    tmpoutfile = tmp2

                else:
                    tmpoutfile = tmp1

                # Make directory if it doesn't exist
                os.makedirs(os.path.dirname(tmpoutfile), exist_ok=True)

                print(command)
                exit_code = os.system(command)
                if exit_code != 0:
                    print(f"ERROR::Command failed with exit code {exit_code}")
                    sys.exit(exit_code)

        # Make directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        if "root://" in output_file:
            print("xrdcp %s %s" % (tmpinfile, output_file))
            os.system("xrdcp %s %s" % (tmpinfile, output_file))
        else:
            print("cp %s %s" % (tmpinfile, output_file))
            os.system("cp %s %s" % (tmpinfile, output_file))

        # Clean the temporary files
        print("rm -f %s" % tmp1)
        os.system("rm -f %s" % tmp1)
        print("rm -f %s" % tmp2)
        os.system("rm -f %s" % tmp2)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()
    PIDCorr(**vars(args))
