'''
Author       : Jie Wu j.wu@cern.ch
Date         : 2024-07-22 05:33:10 +0200
LastEditors  : Jie Wu j.wu@cern.ch
LastEditTime : 2024-07-22 11:12:27 +0200
FilePath     : PIDCorr.py
Description  : 

Copyright (c) 2024 by everyone, All Rights Reserved. 
'''

import os
import sys
import yaml
import argparse


def read_from_yaml(mode, selection_files):
    selection_dict = {}
    for file in selection_files:
        with open(file, 'r') as stream:
            selection_dict.update(yaml.safe_load(stream)[mode])
    return selection_dict


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', help='Path to the input file')
    parser.add_argument('--input-tree-name', help='Name of the tree')
    parser.add_argument('--output-file', help='Output ROOT file')
    parser.add_argument('--data-set', help='Mag and Year, for example MagUp_2016')
    parser.add_argument('--mode', help='Name of the selection in yaml')
    parser.add_argument('--simversion', help='Simversion, where to read samples from')
    parser.add_argument('--tracks-file', nargs='+', help='Yaml file with tracks')
    parser.add_argument('--config-file', help='Path to the PID config sample definition file')
    parser.add_argument('--tmp1', help='Temporary file to apply PID 1st step')
    parser.add_argument('--tmp2', help='Temporary file to apply PID 2nd step')
    return parser


def PIDCorr(input_file, input_tree_name, output_file, data_set, mode, simversion, tracks_file, config_file, tmp1, tmp2):
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

    # IF ON LXPLUS: if /tmp exists and is accessible, use for faster processing
    # IF NOT: use /tmp if you have enough RAM
    # temp_folder = '/tmp'
    # ELSE: use current folder

    ## END OF CONFIG

    # make sure we don't overwrite local files and prefix them with random strings
    import string
    import random

    rand_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))  # get 10 random chars for temp_file prefix

    output_tree = input_tree.split("/")[-1]
    treename = input_tree

    for input_file, output_file, dataset in files:
        tmpinfile = input_file
        tmpoutfile = tmp1
        for track, subst in tracks.items():
            # for var, config in subst.items():
            for var, oldpidvar in subst.items():
                oldpidvar = var if oldpidvar == '' else oldpidvar

                # command = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDCorr.py"
                # command = "python PIDCorr_git.py" #Copied from git, to make local edits
                command = "python scripts/PIDCorr_git_edit.py"  # Copied from git, to make local edits
                command += " -m %s_%s" % (track, ptvar)
                if etavar:
                    command += " -e %s_%s" % (track, etavar)
                elif pvar:
                    command += " -q %s_%s" % (track, pvar)
                else:
                    print('Specify either ETA or P branch name per track')
                    sys.exit(1)
                command += " -n %s" % ntrvar
                command += " -t %s" % treename
                command += " -p %s_%s_PIDCorr" % (track, var)
                command += " -s %s_%s" % (track, oldpidvar)
                command += " -c %s" % config_file
                command += " -d %s" % dataset
                command += " -i %s" % tmpinfile
                command += " -o %s" % tmpoutfile
                command += " -S %s" % simversion

                treename = output_tree
                tmpinfile = tmpoutfile
                if tmpoutfile == tmp1:
                    tmpoutfile = tmp2

                else:
                    tmpoutfile = tmp1

                # Make directory if it doesn't exist
                os.makedirs(os.path.dirname(tmpoutfile), exist_ok=True)

                print(command)
                os.system(command)

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
