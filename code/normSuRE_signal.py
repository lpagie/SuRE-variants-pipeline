import sys
import numpy as np
import pyBigWig
import argparse
import os.path
import glob
import itertools
from os import listdir
import merge_bigwig
import normSuRE_cov_bigwig as normSuRE
#from argparse import Namespace

VERBOSE=True
strands = ''# , '.plus', '.minus'

NORM_SCRIPT="/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP161128_Sure_pipeline_snakemake/code/normSuRE-cov-bigwig.py"

def parse_options():
    # parse user options:

    parser = argparse.ArgumentParser(
        description="generate a normalized SuRE score bigwig for a SuRE library, for plus/minus/combined strands, for a specified set of cDNA replicates")

    parser.add_argument("--sureLib",
                        required=True,
                        help=('Label of SuRE library'))

    parser.add_argument("--sureDataDir",
                        required=True,
                        dest='dataDir',
                        help=('pipeline output directory containing input data (eg LP180619_SuRE42_2/pipelineOutput)'))

    parser.add_argument("--cDNA-groups",
                        required=True,
                        nargs='+',
                        help=('strng(s) which specify the group(s) of cDNA samples (eg, "K562 Hepg2"'))

    parser.add_argument("--useMergedIPCR",
                        action='store_true',
                        help=("Don't merge iPCR bigwigs, use existing merged bigwigs [default: false]"))

    parser.add_argument("--out", 
                        required=True,
                        help=("Filename for output bigwig with normalized SuRE score data")) 

    parser.add_argument("-v", "--verbose", 
                        required=False,
                        action='store_true',
                        help=("print verbose output"))

    options = parser.parse_args()
    return options


def merge_cdna_bigwigs(options):
    dataDir = os.path.join(options.dataDir, 'cDNA')
    for gr in options.cDNA_groups:
        _reps = listdir(dataDir)
        reps = list(itertools.compress(_reps, [gr in s for s in _reps]))
        for rep in reps:
            m_ops=argparse.Namespace(inbw=glob.glob(os.path.join(dataDir, rep, rep+'_chr*.cov.bw')),
                                     outfname=os.path.join(dataDir, rep, rep+'.cov.bw'))
            merge_bigwig.main(m_ops)
            if options.verbose:
                print ('done: merge cDNA bigwigs in '+rep)

            m_ops=argparse.Namespace(inbw=glob.glob(os.path.join(dataDir, rep, rep+'_chr*.cov.plus.bw')),
                                     outfname=os.path.join(dataDir, rep, rep+'.cov.plus.bw'))
            #            merge_bigwig.main(m_ops)
            if options.verbose:
                print ('done: merge plus cDNA bigwigs in '+rep)

            m_ops=argparse.Namespace(inbw=glob.glob(os.path.join(dataDir, rep, rep+'_chr*.cov.minus.bw')),
                                     outfname=os.path.join(dataDir, rep, rep+'.cov.minus.bw'))
            #            merge_bigwig.main(m_ops)
            if options.verbose:
                print ('done: merge minus cDNA bigwigs in '+rep)


def merge_ipcr_bigwigs(options):
    dataDir = os.path.join(options.dataDir, 'iPCR')
    m_ops = argparse.Namespace(inbw=glob.glob(os.path.join(dataDir, 'iPCR_chr*.cov.bw')), 
                               outfname=os.path.join(dataDir, 'iPCR.cov.bw'))
    merge_bigwig.main(m_ops)
    if options.verbose:
        print ('done: merge iPCR bigwigs')

    m_ops=argparse.Namespace(inbw=glob.glob(os.path.join(dataDir, 'iPCR_chr*.cov.plus.bw')),
                             outfname=os.path.join(dataDir, 'iPCR.cov.plus.bw'))
    #        merge_bigwig.main(m_ops)
    if options.verbose:
        print ('done: merge plus iPCR bigwigs')

    m_ops=argparse.Namespace(inbw=glob.glob(os.path.join(dataDir, 'iPCR_chr*.cov.minus.bw')),
                             outfname=os.path.join(dataDir, 'iPCR.cov.minus.bw'))
    #        merge_bigwig.main(m_ops)
    if options.verbose:
        print ('done: merge minus iPCR bigwigs')


def norm_sure_signal(options):
    _reps = listdir(os.path.join(options.dataDir, 'cDNA'))
    for gr in options.cDNA_groups:
        reps = list(itertools.compress(_reps, [gr in s for s in _reps]))
        for strand in strands:
            cDNA_input = [os.path.join(options.dataDir, 'cDNA', r, r+'.cov'+strand+'.bw') for r in reps]
            iPCR_input = os.path.join('iPCR','iPCR.cov'+strand+'.bw')
            outfname = '_'.join([options.sureLib] + reps + ['.cov'+strand+'.bw'])
            print(outfname)
            n_options = argparse.Namespace(out=outfname, cdna=' '.join(cDNA_input), 
                                           input=iPCR_input, pseudo=0)
            print(n_options)
            normSuRE.main(n_options)


def main(options):
    merge_cdna_bigwigs(options)
    if not options.useMergedIPCR: 
        merge_ipcr_bigwigs(options)
    norm_sure_signal(options)
    ## all chromosome names, in correct order
    # CHRS = ['chr'+str(i) for i in (list(range(1,22)) + ['X', 'Y','M'])]
    # CHRS = ['chr'+str(i) for i in (list(range(19,21)))]
#    CHR = []

#    cov = norm_chr(options.input, options.cdna, CHR, options.pseudo)
#    fname = os.path.join(options.out)
#    write_bw(cov, fname, CHR)


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))

    options = parse_options()

    main(options)

