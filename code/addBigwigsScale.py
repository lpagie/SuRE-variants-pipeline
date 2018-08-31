# AUTHOR / DATE
#   Ludo Pagie; 180831, addBigwigsScale.py

# INTRO / BACKGROUND
#   python script to combine two or more bigwigs (intended for normalized SuRE
#   data bigwigs). The data of the bigwigs is summed and then scaled such that
#   the mean of the data is 1.
#   If the bigwigs use different sets of regions this can be specified using
#   option '--libs-differ'). In this case during data import the data is
#   expanded to base resolution. After all data is imported the data again is
#   run length encoded and exported to bigwig.
#
# USAGE / INPUT / ARGUMENTS / OUTPUT
# USAGE:
#   required:
#     --bigwigs: 2 or more bigwig filenames
#     --outfname: bigwig filename for output
#   optional:
#     --chrs: chromosome names to be processed (default: chr1-chr22, chrX)
#     --libs-differ: boolean indicating the bigwig data originates from different libraries
#     --verbose: print debug info
# INPUT:
#   bigwig files (eg normalized SuRE data of different cDNA data sets)
# OUTPUT:
#   bigwig with combined and scaled data

# VERSIONS:
#   

# TODO
#   

import sys
import numpy as np
import pyBigWig
import argparse
import itertools

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
                        default=["chr%s" % c for c in list(range(1,23))+['X']],
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

    return options


def init(bw_fname, chr_lengths):
    # initialize some data dependent variables

    # set chr_lengths from 1st bigwig
    inpt_bw=pyBigWig.open(bw_fname[0])
    for c in chr_lengths:
        chr_lengths[c] = inpt_bw.chroms(c)
    inpt_bw.close()


def read_bigwigs(bw_fname, cov, options):
    # read data from all bigwigs, store the summed data in dictionary cov
    if options.VERBOSE:
        print("start read_bigwigs(...)")
    for f in bw_fname:
        # loop over bigwigs
        inpt_bw=pyBigWig.open(f)
        # loop over chromosomes
        for c in cov:
            inpt_ints = inpt_bw.intervals(c) # access intervals in bigwig
            if cov[c] is None:
                # first bigwig creates complete dictionary element for cov[c]
                cov[c] = {'val':np.array([i[2] for i in inpt_ints]),
                          'start':np.array([i[0] for i in inpt_ints]),
                          'end':np.array([i[1] for i in inpt_ints]),
                          'header':(c, inpt_bw.chroms(c))}
            else:
                # remaining bigwigs only add data to cov[c]['val']
                cov[c]['val'] += np.array([i[2] for i in inpt_ints])
        inpt_bw.close()

    if options.VERBOSE:
        print("exit read_bigwigs(...)")



def read_bigwigs_expand2base(bw_fname, cov, chr_lengths, options):
    # same as previous function (read_bigwigs) but this one deals with bigwigs
    # based on different intervals. All data is expanded to basepair resolution
    # and processed in a chromosome length array.
    if options.VERBOSE:
        print("start read_bigwigs_expand2base(...)")

    for c in cov:
        # loop over chromosomes in outer loop to minimize memory usage
        # per chromosome initialize the chromosome length data vector to 0's,
        # sum data from all bigwigs, run length encode the resulting data
        # vector, and store in a cov dictionary element

        vals = np.zeros(chr_lengths[c]) # init data vector
        for f in bw_fname:
            # loop over bigwigs
            inpt_bw=pyBigWig.open(f)
            inpt_ints = inpt_bw.intervals(c)
            for i in inpt_ints:
                # add data onto base resolution data vector
                vals[i[0]:i[1]] += [i[2]]*(i[1]-i[0])
            inpt_bw.close()
        # run length encode the data vector
        rle = [[k,len(list(i))] for k,i in itertools.groupby(vals)]
        # cumulative sum of rle.lengths for start/end coordinates
        lens_cum = np.cumsum([e[1] for e in rle])
        start = np.concatenate((np.array([0]), lens_cum[:-1]))
        end = lens_cum
        # extract values from rle
        val = np.array([e[0] for e in rle])
        # create dictionary element for 'cov'
        cov[c] = {'val':val, 'start':start, 'end':end, 'header':(c, chr_lengths[c])}

    if options.VERBOSE:
        print("exit read_bigwigs_expand2base(...)")


def scale_coverage(cov):
    # scale the imported data to mean(data)==1
    len_tot = 0. # length of regions in bigwig
    cov_tot = 0. # values of regions in bigwig
    for c in cov:
        # compute total sum of values and length
        ll = cov[c]['end'] - cov[c]['start']
        cov_tot += sum(ll * cov[c]['val'])
        len_tot += sum(ll)

    # compute genome-wide mean score and scale normalized score to mean = 1
    mean = cov_tot/len_tot
    for c in cov:
        cov[c]['val'] = cov[c]['val']/mean

    return cov


def write_bw(cov, fname):
    # cov is a dictionary where each element is a dictionary with 2 lists (start, end),
    # a tuple (header), and 1 np.array (val)
    # the bigwig file needs a header composed of all headers in the dict 'cov'
    # the data is written per chromosome
    # fname is the name of the bigwig output file

    # construct bigwig header from cov
    header = [cov[c]["header"] for c in cov]
    print("writing data to output")
    # open bigwig, write mode
    norm_bw = pyBigWig.open(fname, "w")
    assert(norm_bw is not None)
    # write header to bigwig file
    norm_bw.addHeader(header)
    # loop over elements in 'cov'. for each chromosome in 'cov' write data to bigwig
    for c in cov:
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
    cov = {c:None for c in CHR}
    if options.VERBOSE:
        print('starting init()')
        print(CHR)
        print(options.bigwigs)
        print(chr_lengths)
    init(options.bigwigs, chr_lengths)
    if options.VERBOSE:
        print(chr_lengths)
        print(CHR)
        print(cov)
        print('finished init(..) and starting read_bigwigs(..)')
    if options.libs_differ:
        read_bigwigs_expand2base(options.bigwigs, cov, chr_lengths, options)
    else:
        read_bigwigs(options.bigwigs, cov, options)
    if options.VERBOSE:
        print('finished read_bigwigs(..) and starting scale_coverage(..)')
    scale_coverage(cov)
    if options.VERBOSE:
        print('finished scale_coverage(..) and starting write_bw(..)')
    write_bw(cov, options.outfname)
    if options.VERBOSE:
        print('finished write_bw(..), finished adding and scaling bigwigs')


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))

    options = parse_options()

    main(options)

