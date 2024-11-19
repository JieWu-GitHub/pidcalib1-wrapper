'''
Author       : Jie Wu j.wu@cern.ch
Date         : 2024-11-06 07:33:16 +0100
LastEditors  : Jie Wu j.wu@cern.ch
LastEditTime : 2024-11-06 07:48:44 +0100
FilePath     : PIDGenUtils_edit.py
Description  : 

Copyright (c) 2024 by everyone, All Rights Reserved. 
'''

###############################################################################
# (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
'''Common utils for PIDGen and PIDCorr.'''
import argparse
from ROOT import gROOT, TFile, TTree

gROOT.ProcessLine(
    "struct MyStruct { \
Double_t newpid; \
Double_t hint; \
Double_t hintmc; \
}"
)

# Defaults for commandline arguments
defaults = {
    'tree': 'tree',
    'output': 'output.root',
    'pidvar': 'PID_gen',
    'ptvar': 'Pt',
    'pvar': 'P',
    'etavar': None,
    'ntrvar': 'nTracks',
    'lowerpid': None,
    'config': 'p_V3ProbNNp',
    'dataset': 'MagDown_2011',
    'var': 'default',
    'seed': None,
    'ntrscale': None,
    'calibstat': False,
    'noclone': False,
    'simpidvar': 'PID',
    'simversion': 'sim08',
    'outtree': None,
    'localrootdir': None,
}


def get_argparser(pidgen=True):
    '''Get the ArgumentParser for PIDGen or PIDCorr.'''
    parser = argparse.ArgumentParser(description='PIDGen' if pidgen else 'PIDCorr')
    parser.add_argument('-i', '--input', type=str, default=None, help='Input file name')
    parser.add_argument('-t', '--tree', type=str, default=defaults['tree'], help='Input tree name')
    parser.add_argument('-o', '--output', type=str, default=defaults['output'], help='Output file name')
    parser.add_argument('-p', '--pidvar', type=str, default=defaults['pidvar'], help='Output name for the generated/corrected PID variable')
    parser.add_argument('-m', '--ptvar', type=str, default=defaults['ptvar'], help='Pt variable')
    parser.add_argument('-q', '--pvar', type=str, default=defaults['pvar'], help='P variable')
    parser.add_argument('-e', '--etavar', type=str, default=defaults['etavar'], help='Eta variable (if None, calculated from P and Pt)')
    parser.add_argument('-n', '--ntrvar', type=str, default=defaults['ntrvar'], help='Ntracks variable')
    parser.add_argument('-l', '--lowerpid', type=str, default=defaults['lowerpid'], help='Lower PID value to generate')
    parser.add_argument('-c', '--config', type=str, default=defaults['config'], help=('PID response to sample. Run without giving -i or --input to' ' list available configs.'))
    parser.add_argument('-d', '--dataset', type=str, default=defaults['dataset'], help='Dataset (polarity_year)')
    parser.add_argument('-v', '--var', type=str, default=defaults['var'], help='Variation (default, syst_N, stat_N etc.)')
    parser.add_argument('-f', '--ntrscale', type=float, default=defaults['ntrscale'], help='Scale factor for nTracks variable (default - no scaling)')
    parser.add_argument('-a', '--calibstat', help='Add calibration statistics branch', action="store_const", const=True, default=defaults['calibstat'])
    parser.add_argument('--noclone', help='Don\'t clone the original tree in the output.', action='store_true')
    parser.add_argument('--outtree', help='Name of the output TTree', default=defaults['outtree'])
    parser.add_argument(
        '--localrootdir',
        help='Use local directory with input source calibration ROOT files, instead of the centrally stored on in the eos. Should be down to the folder /lhcb/wg/PID/PIDGen',
        default=defaults['localrootdir'],
    )

    if pidgen:
        parser.add_argument('-s', '--seed', type=str, default=defaults['seed'], help='Initial random seed')
    else:
        parser.add_argument(
            '-s', '--simpidvar', type=str, default=defaults['simpidvar'], help=('Original, simulated PID variable to correct. ' 'Eg, <head>_PID(K|p|mu|e) or <head>_ProbNN(pi|k|p|mu|e)')
        )
        parser.add_argument('-S', '--simversion', type=str, default=defaults['simversion'], help='Simulation version ("sim08" or "sim09" for Run1, "run2" for Run2)')

    return parser


def make_output_tree(tree, noclone, outfilename, outtree):
    '''Make the output TTree.'''
    if isinstance(outtree, TTree):
        newtree = outtree
        outfile = None
    else:
        outfile = TFile.Open(outfilename, "RECREATE")
        if noclone:
            newtree = TTree(tree.GetName(), tree.GetTitle())
        else:
            newtree = tree.CloneTree(0)
        if outtree:
            newtree.SetName(outtree)
    return newtree, outfile


def get_fill_objects(newtree, outtree, branches):
    '''Get objects to call Fill on.'''
    if isinstance(outtree, TTree):
        return branches
    return [newtree]
