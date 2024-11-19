'''
Author       : Martin Andersson and Jie Wu j.wu@cern.ch and 
Date         : 2024-07-20 16:21:32 +0200
LastEditors  : Jie Wu j.wu@cern.ch
LastEditTime : 2024-11-07 07:28:37 +0100
FilePath     : PIDCorr_edit.py
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
from __future__ import print_function
from __future__ import division
from past.utils import old_div
import os, sys, math
from ROOT import TFile, TH1F, OneDimPhaseSpace, CombinedPhaseSpace, BinnedDensity, Logger, TTree, MyStruct, addressof, std, gRandom
import ROOT as r
from math import sqrt, log
from dataclasses import dataclass

import yaml


# from PIDPerfScripts.PIDGenUtils import get_argparser, defaults, make_output_tree, get_fill_objects
import PIDGenExpert.Run1.Config as ConfigRun1
import PIDGenExpert.Run2.Config as ConfigRun2
import PIDGenExpert.Run1.ConfigMC as ConfigMCSim08
import PIDGenExpert.Run1.ConfigMCSim09 as ConfigMCSim09
import PIDGenExpert.Run2.ConfigMC as ConfigMCRun2


from PIDGenUtils_edit import get_argparser, defaults, make_output_tree, get_fill_objects


def read_from_yaml(selection_files, mode=None):
    selection_dict = {}
    selection_files = selection_files.split(';') if isinstance(selection_files, str) else selection_files

    for file in selection_files:
        with open(file, 'r') as stream:
            selection_dict.update(yaml.safe_load(stream)[mode]) if mode else selection_dict.update(yaml.safe_load(stream))
    return selection_dict


def check_file_existence(url: str) -> bool:
    """
    Check if a ROOT file exists and can be opened.

    Args:
    url (str): URL or path to the ROOT file.

    Returns:
    bool: True if the file exists and can be opened, False otherwise.
    """

    if not url:
        return False

    try:

        # Suppress ROOT error messages temporarily
        old_error_level = r.gErrorIgnoreLevel
        r.gErrorIgnoreLevel = r.kFatal

        # Attempt to open the file
        input_file = TFile.Open(url, "READ")

        # Restore original error level
        r.gErrorIgnoreLevel = old_error_level

        # Check if the file is opened successfully and is not a zombie
        if input_file and not input_file.IsZombie():
            input_file.Close()
            return True
        else:
            return False

        # # Attempt to open the file
        # input_file = TFile.Open(url, "READ")

        # # Check if the file is opened successfully and is not a zombie
        # if input_file and not input_file.IsZombie():
        #     input_file.Close()
        #     return True
        # else:
        #     return False

    except Exception as e:
        # Optionally print out the error for debugging
        # print(f"Error opening file: {e}")
        return False


# confdict links TRUEID with the correct sample
# ------------------------------------
confdict = {
    '13': 'mu',
    '211': 'pi',
    '321': 'K',
    '2212': 'p',
    '11': 'e',
}


@dataclass
class COUNTERS:
    # From PID maps
    counter_nocalib: int = 0
    counter_nomc: int = 0
    counter_lowcalib: int = 0
    counter_lowmc: int = 0

    # From the tree
    n_mu: int = 0
    n_pi: int = 0
    n_k: int = 0
    n_p: int = 0
    n_e: int = 0
    n_branch: int = 0

    n_mu_nocalib: int = 0
    n_pi_nocalib: int = 0
    n_k_nocalib: int = 0
    n_p_nocalib: int = 0
    n_e_nocalib: int = 0
    n_branch_nocalib: int = 0

    count: int = 0
    nPIDGen: int = 0
    nPIDCorr: int = 0


# ------------------------------------


def main():
    '''Main function to run PIDCorr.'''
    parser = get_argparser(False)
    parser.print_help()
    args = parser.parse_args()

    print("PIDCorr arguments:\n", args)

    infilename = args.input
    intree = args.tree
    outfilename = args.output
    pidvar = args.pidvar
    ptvar = args.ptvar
    pvar = args.pvar
    etavar = args.etavar
    ntrvar = args.ntrvar
    minpid = args.lowerpid
    oldpidvar = args.simpidvar
    config_file = args.config
    dataset = args.dataset
    variant = args.var
    simversion = args.simversion
    outtree = args.outtree
    ntrscale = args.ntrscale
    addcalibstat = args.calibstat
    localrootdir = args.localrootdir

    # The name of the track
    track = oldpidvar.split('_')[0]

    # the config sample name
    PIDconfKey = pidvar.split('_')[1]

    if not infilename:
        print("Usage: PIDCorr.py [options]")
        #  print "  For the usage example, look at pid_transform.sh file"
        print("  Available PID configs for Run1/sim08 are: ")
        for i in sorted(ConfigRun1.configs.keys()):
            if i in list(ConfigMCSim08.configs.keys()):
                print("    ", i)
        print("  Available PID configs for Run1/sim09 are: ")
        for i in sorted(ConfigRun1.configs.keys()):
            if i in list(ConfigMCSim09.configs.keys()):
                print("    ", i)
        print("  Available PID configs for Run2/sim09 are: ")
        for i in sorted(ConfigRun2.configs().keys()):
            if i in list(ConfigMCRun2.configs.keys()):
                print("    ", i)

        # Exit politely
        sys.exit(0)

    # Read in the config file
    PIDconfig = read_from_yaml(config_file)
    calibOption = PIDconfig['option']
    calibStrict = PIDconfig['strict']
    samples = PIDconfig['samples']

    run_pid_corr(
        infilename=infilename,
        intree=intree,
        outfilename=outfilename,
        pidvar=pidvar,
        ptvar=ptvar,
        pvar=pvar,
        etavar=etavar,
        ntrvar=ntrvar,
        minpid=minpid,
        oldpidvar=oldpidvar,
        # config=conf,
        samples=samples,
        calibOption=calibOption,
        calibStrict=calibStrict,
        dataset=dataset,
        variant=variant,
        simversion=simversion,
        ntrscale=ntrscale,
        addcalibstat=addcalibstat,
        noclone=args.noclone,
        outtree=outtree,
        PIDconfKey=PIDconfKey,
        track=track,
        localrootdir=localrootdir,
    )


def run_pid_corr(
    infilename=None,
    intree=defaults['tree'],
    outfilename=defaults['output'],
    pidvar=defaults['pidvar'],
    ptvar=defaults['ptvar'],
    pvar=defaults['pvar'],
    etavar=defaults['etavar'],
    ntrvar=defaults['ntrvar'],
    minpid=defaults['lowerpid'],
    oldpidvar=defaults['simpidvar'],
    # config=defaults['config'],
    samples=None,
    calibOption=None,
    calibStrict=True,
    dataset=defaults['dataset'],
    variant=defaults['var'],
    simversion=defaults['simversion'],
    ntrscale=defaults['ntrscale'],
    addcalibstat=defaults['calibstat'],
    noclone=defaults['noclone'],
    outtree=defaults['outtree'],
    PIDconfKey="None",  # Added
    track="None",
    localrootdir=defaults['localrootdir'],
):  # Added
    '''Run PIDCorr with the given config. intree can be a TTree or the name of
    the tree to be retrieved from the file named infilename. Similarly, outtree
    can either be the name of the tree to be saved to the file outfilename
    (default: the same name as intree), or a TTree to which the PIDGen branches
    will be added.'''
    if not infilename and not isinstance(intree, TTree):
        raise ValueError("If intree isn't a TTree then you must specify infilename!")

    if variant == "default":
        variant = "distrib"  # to do: change this name in CreatePIDPdf

    # datapdf = Config.eosrootdir + "/" + config + "/" + dataset + "_" + variant + ".root"
    # simpdf = ConfigMC.eosrootdir + "/" + config + "/" + dataset + "_" + variant + ".root"

    # Set year and run
    year = None
    run = None
    try:
        year = dataset.split("_")[1]
    except:
        print('Dataset format "%s" not recognized. Should be {MagUp,MagDown}_[Year]' % dataset)
        quit()
    if year in ["2011", "2012"]:
        run = 'run1'
    elif year in ["2015", "2016", "2017", "2018"]:
        run = 'run2'
    else:
        print('Data taking year "%s" not recognized' % year)
        quit()

    # Set the correct Config and ConfigMC (for PIDCorr only)
    if calibOption == "PIDCorr":
        if simversion == "sim08":
            ConfigMC = ConfigMCSim08
            Config = ConfigRun1
        elif simversion == "sim09":
            ConfigMC = ConfigMCSim09
            Config = ConfigRun1
        elif simversion == "run2":
            ConfigMC = ConfigMCRun2
            Config = ConfigRun2
        else:
            print("Simulation version %s unknown" % simversion)
            sys.exit(1)
    elif calibOption == "PIDGen":
        if run == 'run1':
            Config = ConfigRun1
        elif run == 'run2':
            Config = ConfigRun2

    # Read the input tree
    if isinstance(intree, TTree):
        tree = intree
        infile = None
    else:
        infile = TFile.Open(infilename)
        tree = infile.Get(intree)
    if not tree:
        print("Ntuple not found!")
        sys.exit(1)

    nentries = tree.GetEntries()

    # datakde = BinnedDensity("KDEPDF", phsp, datapdf)
    # simkde = BinnedDensity("KDEPDF", phsp, simpdf)

    s = MyStruct()

    newtree, outfile = make_output_tree(tree, noclone, outfilename, outtree)

    branches = [newtree.Branch(pidvar, addressof(s, "newpid"), f'{pidvar}/D')]
    addcalibstat = True
    if addcalibstat:
        branches.append(newtree.Branch(f"{pidvar}_calibstat", addressof(s, "hint"), f'{pidvar}_calibstat/D'))
        branches.append(newtree.Branch(f"{pidvar}_mcstat", addressof(s, "hintmc"), f'{pidvar}_mcstat/D')) if calibOption == "PIDCorr" else None
    fillobjs = get_fill_objects(newtree, outtree, branches)

    if infile:
        infile.cd()

    # Just setting some range for the TH1F as this will change later
    # ----------------------------------
    minpid = -150
    pidmax = 150
    # ----------------------------------
    hdata = TH1F("hdata", "h", 100, minpid, pidmax)
    hsim = TH1F("hsim", "h", 100, minpid, pidmax)

    # print(transform_backward)

    if noclone:
        tree.SetBranchStatus('*', False)
        tree.SetBranchStatus(oldpidvar, True)
        tree.SetBranchStatus(ptvar, True)
        tree.SetBranchStatus(ntrvar, True)
    var_code = compile("i.%s" % oldpidvar, '<string>', 'eval')
    # pid_code = compile(transform_backward, '<string>', 'eval')
    # oldpid_code = compile(transform_forward, '<string>', 'eval')
    log_pt_code = compile("log(i.%s)" % ptvar, '<string>', 'eval')
    if etavar == None:
        p_code = compile("i.%s" % pvar, '<string>', 'eval')
        pt_code = compile("i.%s" % ptvar, '<string>', 'eval')
        if noclone:
            tree.SetBranchStatus(pvar, True)
    else:
        eta_code = compile("i.%s" % etavar, '<string>', 'eval')
        if noclone:
            tree.SetBranchStatus(etavar, True)
    if ntrscale:
        log_ntracks_code = compile("log(float(i.%s)*%f)" % (ntrvar, ntrscale), '<string>', 'eval')
    else:
        log_ntracks_code = compile("log(float(i.%s))" % ntrvar, '<string>', 'eval')

    Logger.setLogLevel(1)

    # Load all samples
    # -------------------------
    datakdes = {}
    simkdes = {}

    # ----- for datapdf -----
    for particle_id, particle_name in confdict.items():

        # 1) Check if the particle_PIDconfKey exists in the config file
        particle_PIDconfKey = f'{particle_name}_{PIDconfKey}'
        if particle_PIDconfKey not in samples[run]:
            print(f"INFO: Combination of {particle_name} and {PIDconfKey} is not corrected or resampled because the key do not exist in the configuration file.\n")
            continue
        # 2) Check if the datapdf exists
        datapdf = Config.eosrootdir + "/" + samples[run][particle_PIDconfKey] + "/" + dataset + "_" + variant + ".root"
        if os.path.exists(localrootdir):  # Check if the local directory exists, if not use eosrootdir
            _datapdf_path_list = datapdf.split('/lhcb/wg/PID/PIDGen')
            datapdf = localrootdir + _datapdf_path_list[-1]

        if check_file_existence(datapdf):
            __configs, phsp, __minpid, __pidmin, __pidmax, __transform_forward, __transform_backward = get_settings_from_sample(samples[run][f'{particle_name}_{PIDconfKey}'], run, dataset, variant)
            datakdes[particle_name] = BinnedDensity("KDEPDF", phsp, datapdf)
            print(f"INFO: Set datakde for particle: {particle_name}, path to pdf: {datapdf}")
        else:
            print(f"INFO: Combination of {particle_name} and {PIDconfKey} is not resampled because the file to the pdf do not exists: {datapdf}\n")

    # Check datakdes is not empty
    if not datakdes:
        print("ERROR: No pdfs for data is loaded, exiting...")
        exit(1)

    # ----- for simpdf -----
    if calibOption == "PIDCorr":
        for particle_name in datakdes:
            if calibOption == "PIDCorr":

                particle_PIDconfKey = f'{particle_name}_{PIDconfKey}'

                # Check if the simpdf exists
                simpdf = ConfigMC.eosrootdir + "/" + samples[run][particle_PIDconfKey] + "/" + dataset + "_" + variant + ".root"

                if os.path.exists(localrootdir):  # Check if the local directory exists, if not use eosrootdir
                    _simpdf_path_list = simpdf.split('/lhcb/wg/PID/PIDGen')
                    simpdf = localrootdir + _simpdf_path_list[-1]

                if check_file_existence(simpdf):
                    __configs, phsp, __minpid, __pidmin, __pidmax, __transform_forward, __transform_backward = get_settings_from_sample(
                        samples[run][f'{particle_name}_{PIDconfKey}'], run, dataset, variant
                    )
                    simkdes[particle_name] = BinnedDensity("KDEPDF", phsp, simpdf)
                    print(f"INFO: Set simkde for particle: {particle_name}, path to pdf: {simpdf}")
                else:
                    if calibStrict:
                        print(f'ERROR: Combination of {particle_name} and {PIDconfKey} is not corrected or resampled because the file to the pdf do not exists: {simpdf}')
                        exit(1)
                    else:
                        print(f"INFO: No simulation, so only set datakde for particle (do resampling instead): {particle_name}")

        # Check simkdes is not empty
        if not simkdes:
            print("ERROR: PIDCorr is called but no pdfs for simulation is loaded. Please check the MC pdf samples in the configuration file, or use PIDGen instead.")
            exit(1)

    # -------------------------
    counter = COUNTERS()

    # Set trueid compiler and counters
    # -------------------------
    trueid = compile("i.%s" % track + "_TRUEID", '<string>', 'eval')

    # -------------------------
    for i in tree:
        if calibOption == "PIDCorr":
            resample_corr = True
            resample_gen = False
        else:
            resample_corr = False
            resample_gen = True

        point = std.vector('double')(4)
        # point[0] = (pidmin + pidmax) / 2.
        point[1] = eval(log_pt_code)
        point[3] = eval(log_ntracks_code)
        if etavar == None:
            point[2] = -math.log(math.tan(math.asin(old_div(eval(pt_code), eval(p_code))) / 2.0))
        else:
            point[2] = eval(eta_code)
        # -------------------------
        TRUEID = eval(trueid)
        true_particle = confdict.get(str(abs(TRUEID)), 'unknown')

        if true_particle == "mu":
            counter.n_mu += 1
        elif true_particle == "pi":
            counter.n_pi += 1
        elif true_particle == "K":
            counter.n_k += 1
        elif true_particle == "p":
            counter.n_p += 1
        elif true_particle == "e":
            counter.n_e += 1
        else:
            counter.n_branch += 1

        # Check if the particle_PIDconfKey exists in the config file
        if true_particle in datakdes:

            sample = samples[run][f'{true_particle}_{PIDconfKey}']

            datakde = datakdes[true_particle]
            if resample_corr:
                if true_particle in simkdes:
                    simkde = simkdes[true_particle]
                else:
                    resample_corr = False
                    resample_gen = True

            __configs, __phsp, minpid, pidmin, pidmax, transform_forward, transform_backward = get_settings_from_sample(sample, run, dataset, variant)
            # print("-----------------------------------------------------------")
            # print("-----------------------------------------------------------")
            # print("DEBUG: Settings for sample", sample)
            # print("DEBUG: run", run)
            # print("DEBUG: dataset", dataset)
            # print("DEBUG: variant", variant)
            # print("DEBUG: pidmax", pidmax)
            # print("DEBUG: pidmin", pidmin)
            # print("DEBUG: minpid", minpid)
            # print("-----------------------------------------------------------")
            # print("-----------------------------------------------------------")
            pid_code = compile(transform_backward, '<string>', 'eval')
            oldpid_code = compile(transform_forward, '<string>', 'eval')
            hdata.SetAxisRange(minpid, pidmax)
            hdata.SetBins(100, minpid, pidmax)
            hsim.SetAxisRange(minpid, pidmax)
            hsim.SetBins(100, minpid, pidmax)

            point[0] = (pidmin + pidmax) / 2.0
            # -------------------------

            # print(point[0], point[1], point[2], point[3])
            # exit(1)

            hdata.Reset()
            hsim.Reset()

            # datakde, common for both resampling methods
            datakde.slice(point, 0, hdata)
            s.hint = hdata.Integral()
            if s.hint == 0:
                counter.counter_nocalib += 1
            elif s.hint < 10:
                counter.counter_lowcalib += 1

            # Select the correct resampling method
            if resample_corr:

                # simkde, only for PIDCorr
                simkde.slice(point, 0, hsim)
                s.hintmc = hsim.Integral()
                if s.hintmc == 0:
                    counter.counter_nomc += 1
                elif s.hintmc < 10:
                    counter.counter_lowmc += 1

                # PIDCorr algorithm
                x = eval(var_code)
                oldpid = x
                if transform_forward == "x" or x >= 0:
                    oldpid = eval(oldpid_code)
                    if oldpid < pidmin or oldpid > pidmax:
                        x = oldpid
                    else:
                        x = datakde.transform(hsim, hdata, oldpid)
                    s.newpid = eval(pid_code)
                else:  # The case for ProbNN<0, just leave as it is
                    s.newpid = x

                # Fill the tree
                for obj in fillobjs:
                    obj.Fill()

                if counter.count % 1000 == 0:
                    print(
                        "Event %d/%d : Pt=%f, Eta=%f, Ntr=%f, OldPID=%f, PIDCorr=%f, X=%f, CalibStat=%f, MCStat=%f"
                        % (counter.count, nentries, point[1], point[2], point[3], oldpid, s.newpid, x, s.hint, s.hintmc)
                    )

                counter.nPIDCorr += 1

            elif resample_gen:

                # PIDGen algorithm
                x = eval(var_code)
                oldpid = x

                if hdata.Integral() > 0:
                    x = hdata.GetRandom()
                else:
                    x = minpid + (pidmax - minpid) * gRandom.Rndm()

                s.newpid = eval(pid_code)

                # Fill the tree
                for obj in fillobjs:
                    obj.Fill()

                if counter.count % 1000 == 0:  # or dontResamp) :
                    print("PIDGen Event %d/%d : Pt=%f, Eta=%f, Ntr=%f, OldPID=%f, PIDGen=%f, X=%f" % (counter.count, nentries, point[1], point[2], point[3], oldpid, s.newpid, x))

                counter.nPIDGen += 1

        else:
            if true_particle == "mu":
                counter.n_mu_nocalib += 1
            elif true_particle == "pi":
                counter.n_pi_nocalib += 1
            elif true_particle == "K":
                counter.n_k_nocalib += 1
            elif true_particle == "p":
                counter.n_p_nocalib += 1
            elif true_particle == "e":
                counter.n_e_nocalib += 1
            else:
                counter.n_branch_nocalib += 1

            s.newpid, s.hint, s.hintmc = -9999.0, -9999.0, -9999.0

            # Fill the tree
            for obj in fillobjs:
                obj.Fill()

            if counter.count % 1000 == 0:
                print("Event %d/%d : Pt=%f, Eta=%f, Ntr=%f, true_particle=%s" % (counter.count, nentries, point[1], point[2], point[3], true_particle))

            # print(f"INFO: Particle {true_particle} is not corrected or resampled because the key do not exist in the configuration file.\n")

        counter.count += 1

    if noclone:
        tree.SetBranchStatus('*', True)

    # If adding branches to an existing output tree, make sure
    # the number of entries is set correctly for the tree.
    if isinstance(outtree, TTree):
        newtree.SetEntries(tree.GetEntries())

    if infile:
        infile.Close()
    if outfile:
        outfile.cd()
        newtree.Write()
        outfile.Close()

    # Print statistics
    print("PID Resampling finished.")
    print(f"Calibration option:    {calibOption}")
    print(f"    Total number of events in the tree:               {nentries}")
    print(f"    Total number of events processed:                 {counter.count}")
    print(f"Resampled {counter.nPIDCorr} with Correlations (PID transformation) and {counter.nPIDGen} without (PIDGen)")

    # Calculating values for cleaner code
    mu_calib = counter.n_mu - counter.n_mu_nocalib
    pi_calib = counter.n_pi - counter.n_pi_nocalib
    k_calib = counter.n_k - counter.n_k_nocalib
    p_calib = counter.n_p - counter.n_p_nocalib
    e_calib = counter.n_e - counter.n_e_nocalib
    n_branch_calib = counter.n_branch - counter.n_branch_nocalib

    # Creating a consistent format for items with slash
    format_str = "    {:<40} {:>10d} / {:<10d}"
    print(format_str.format("Total number of calibrated muons:", mu_calib, counter.n_mu))
    print(format_str.format("Total number of calibrated pions:", pi_calib, counter.n_pi))
    print(format_str.format("Total number of calibrated kaons:", k_calib, counter.n_k))
    print(format_str.format("Total number of calibrated protons:", p_calib, counter.n_p))
    print(format_str.format("Total number of calibrated electrons:", e_calib, counter.n_e))
    print(format_str.format("Total number of events from branch:", n_branch_calib, counter.n_branch))

    print("Events calib stats:")
    print(f"    Events with no calibration (n==0):                {counter.counter_nocalib}")
    print(f"    Events with low calib. stats (0<n<10):            {counter.counter_lowcalib}")
    print(f"    Events with no MC (n==0):                         {counter.counter_nomc}")
    print(f"    Events with low MC stats (0<n<10):                {counter.counter_lowmc}")


def get_settings_from_sample(sample, run, dataset, variant, minpid=None):
    if run == 'run1':
        # print("INFO: Getting Run1 settings")
        # calibfilename = ConfigRun1.eosrootdir + "/" + sample + "/" + "%s_%s.root" % (dataset, variant)
        transform_forward = ConfigRun1.configs[sample]['transform_forward']
        transform_backward = ConfigRun1.configs[sample]['transform_backward']
        configs = ConfigRun1.configs
    elif run == 'run2':
        # print("INFO: Getting Run2 settings")
        # calibfilename = ConfigRun2.eosrootdir + "/" + sample + "/" + "%s_%s.root" % (dataset, variant)
        configs = ConfigRun2.configs()
        if 'gamma' in list(configs[sample].keys()):
            gamma = configs[sample]['gamma']
            if gamma < 0:
                transform_forward = "(1.-(1.-x)**%f)" % abs(gamma)
                transform_backward = "(1.-(1.-x)**%f)" % (1.0 / abs(gamma))
            elif gamma == 1.0:
                transform_forward = "x"
                transform_backward = "x"
            else:
                transform_forward = "((x)**%f)" % abs(gamma)
                transform_backward = "((x)**%f)" % (1.0 / abs(gamma))
        else:
            transform_forward = configs[sample]['transform_forward']
            transform_backward = configs[sample]['transform_backward']
    else:
        print("ERROR: Incorrect run=", run, "EXITING...")
        exit()

    pidmin = 0.0
    pidmax = 1.0
    if 'limits' in configs[sample]:
        pidmin = configs[sample]['limits'][0]
        pidmax = configs[sample]['limits'][1]
    if minpid == None:
        # print("DEBUG: minpid is none")
        minpid = pidmin
    else:
        # print("DEBUG: minpid is not none")
        minpid = float(minpid)
        if minpid < pidmin:
            minpid = pidmin

    # Calculate the minimum PID variable to generate (after transformation)
    x = pidmin
    pidmin = eval(transform_forward)
    x = pidmax
    pidmax = eval(transform_forward)
    x = minpid
    minpid = eval(transform_forward)

    pid_phsp = OneDimPhaseSpace("PIDPhsp", pidmin, pidmax)
    mom_phsp = OneDimPhaseSpace("MomPhsp", 5.5, 9.5)
    eta_phsp = OneDimPhaseSpace("EtaPhsp", 1.5, 5.5)
    ntr_phsp = OneDimPhaseSpace("NtrPhsp", 3.0, 6.5)
    pidmom_phsp = CombinedPhaseSpace("PIDMomPhsp", pid_phsp, mom_phsp)
    pidmometa_phsp = CombinedPhaseSpace("PIDMomEtaPhsp", pidmom_phsp, eta_phsp)
    phsp = CombinedPhaseSpace("FullPhsp", pidmometa_phsp, ntr_phsp)

    return configs, phsp, minpid, pidmin, pidmax, transform_forward, transform_backward


if __name__ == '__main__':
    main()
