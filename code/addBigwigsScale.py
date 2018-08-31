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
        description="Combine bigwigs by adding values and scaling to mean(...)=1")

    parser.add_argument("--bigwigs",
                        nargs='+',
                        required=True,
                        help=('bigwig filenames'))

    parser.add_argument("--outfname", 
                        required=True,
                        help=("Filename for output bigwig with normalized data"))

    parser.add_argument("--chrs",
                        required=False,
                        nargs='+',
                        default=None,
                        help=("Chromosome names to be included in output"))

    parser.add_argument("--libs-differ",
                        required=False,
                        action='store_true',
                        help=("input bigwigs based on different regions; unify"))

    parser.add_argument("-v", "--VERBOSE", 
                        required=False,
                        action='store_true',
                        help=("print verbose output"))

    options = parser.parse_args()

    if options.chrs is None:
        options.chrs = ["chr%s" % c for c in list(range(1,23))+['X']]

    return options


def init(bw_fname, chroms, cov, chr_lengths):
    # import chroms from first input file if chroms is an empty list
    if not chroms:
        inpt_bw=pyBigWig.open(bw_fname[0])
        chroms = list(inpt_bw.chroms().keys())
        chr_lengths = inpt_bw.chroms()
        inpt_bw.close()
        for c in chroms:
            chroms.append(c)

    inpt_bw=pyBigWig.open(bw_fname[0])
    for c in chroms:
        chr_lengths[c] = inpt_bw.chroms(c)
    inpt_bw.close()

    for c in chroms:
        cov[c] = None


def read_bigwigs(bw_fname, CHR, cov, options):
    if options.VERBOSE:
        print("start read_bigwigs(...)")
    for f in bw_fname:
        inpt_bw=pyBigWig.open(f)
#        chrs_in = inpt_bw.chroms().keys()
#        assert (set(chrs_in) == set(CHR)), \
#            "set of chromosomes in bigwig (%s; %s)\ndiffers from used set (%s)" % \
#                    (f, " ".join(chrs_in), " ".join(CHR)) 
        for c in CHR:
            inpt_ints = inpt_bw.intervals(c)
            if cov[c] is None:
                cov[c] = {'val':np.array([i[2] for i in inpt_ints]),
                          'start':np.array([i[0] for i in inpt_ints]),
                          'end':np.array([i[1] for i in inpt_ints]),
                          'header':(c, inpt_bw.chroms(c))}
            else:
                cov[c]['val'] += np.array([i[2] for i in inpt_ints])
        inpt_bw.close()

    if options.VERBOSE:
        print("exit read_bigwigs(...)")



def read_bigwigs_expand2base(bw_fname, CHR, cov, chr_lengths, options):
    if options.VERBOSE:
        print("start read_bigwigs_expand2base(...)")

    for c in CHR:
        inpt_bw=pyBigWig.open(bw_fname[0])
        inpt_ints = inpt_bw.intervals(c)
        vals = np.zeros(chr_lengths[c])
        for f in bw_fname:
            inpt_bw=pyBigWig.open(f)
            inpt_ints = inpt_bw.intervals(c)
            for i in inpt_ints:
                vals[i[0]:i[1]] += [i[2]]*(i[1]-i[0])
            inpt_bw.close()
        rle = [[k,len(list(i))] for k,i in itertools.groupby(vals)]
        lens_cum = np.cumsum([e[1] for e in rle])
        start = np.concatenate((np.array([0]), lens_cum[:-1]))
        end = lens_cum
        val = np.array([e[0] for e in rle])
        cov[c] = {'val':val, 'start':start, 'end':end, 'header':(c, chr_lengths[c])}
        print(start[-10:])
        print(end[-10:])
        print(val[-10:])
        print(c)
        print(chr_lengths[c])
#        print(list(rle.values())[-10:])
#        print(list(rle.keys())[-10:])

    if options.VERBOSE:
        print("exit read_bigwigs_expand2base(...)")


def scale_coverage(cov, CHR):
    len_tot = 0.
    cov_tot = 0.
    for c in CHR:
        ll = cov[c]['end'] - cov[c]['start']
        cov_tot += sum(ll * cov[c]['val'])
        len_tot += sum(ll)

    # compute genome-wide mean score and scale normalized score to mean = 1
    mean = cov_tot/len_tot
    for c in CHR:
        cov[c]['val'] = cov[c]['val']/mean

    return cov


def write_bw(cov, fname, CHR):
    # cov is a dictionary where each element is a dictionary with 2 lists (start, end),
    # a tuple (header), and 1 np.array (norm)
    # the bigwig file needs a header composed of all headers in the dict 'cov'
    # the data is written per chromosome
    # fname is the name of the bigwig output file

    # construct bigwig header from cov
    header = [cov[c]["header"] for c in CHR]
    print("writing data to output")
    # open bigwig, write mode
    norm_bw = pyBigWig.open(fname, "w")
    assert(norm_bw is not None)
    # write header to bigwig file
    norm_bw.addHeader(header)
    # loop over elements in 'cov'. for each chromosome in 'cov' write data to bigwig
    for c in CHR:
        print("exporting data for chrom '%s'" % str(c))
        chrs = [c]*len(cov[c]["start"])
        norm_bw.addEntries(chrs, 
                           cov[c]["start"].tolist(),
                           ends=cov[c]["end"].tolist(),
                           values=cov[c]["val"].tolist())
    # close bigwig file
    print("closing bw")
    norm_bw.close()
    print("writing to output finished")


def main(options):
    CHR = options.chrs
    chr_lengths = {c:None for c in CHR}
    cov = {}
    if options.VERBOSE:
        print('starting init()')
        print(CHR)
        print(options.bigwigs)
        print(chr_lengths)
    init(options.bigwigs, CHR, cov, chr_lengths)
    if options.VERBOSE:
        print(chr_lengths)
        print(CHR)
        print(cov)
        print('finished init(..) and starting read_bigwigs(..)')
    if options.libs_differ:
        read_bigwigs_expand2base(options.bigwigs, CHR, cov, chr_lengths, options)
    else:
        read_bigwigs(options.bigwigs, CHR, cov, options)
    if options.VERBOSE:
        print('finished read_bigwigs(..) and starting scale_coverage(..)')
    scale_coverage(cov, CHR)
    if options.VERBOSE:
        print('finished scale_coverage(..) and starting write_bw(..)')
    write_bw(cov, options.outfname, CHR)
    if options.VERBOSE:
        print('finished write_bw(..), finished adding and scaling bigwigs')


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))

    options = parse_options()

    main(options)

