import sys
import numpy as np
import pyBigWig
import argparse
import os.path
import itertools

def parse_options():
    # parse user options:

    parser = argparse.ArgumentParser(
        description="generate a normalized SuRE score bigwig")
    parser.add_argument("--pseudo", 
                        default=0,
                        help=("Boolean; adding a pseudocount (+1) to input data"))

    parser.add_argument("--out", 
                        required=True,
                        help=("Filename for output bigwig with normalized SuRE score data")) 

    parser.add_argument("--cdna", 
                        required=True,
                        nargs='+', 
                        help=("cDNA input files, ie bigwig files with coverage data"
                              "Multiple names can be given"))
    parser.add_argument("--input", 
                        required=True,
                        nargs='+', 
                        help=("Filenames for input coverage bigwig files (either iPCR or cDNA)")
                       )
    options = parser.parse_args()
    return options

def read_inpt(inpt_fname, CHR, ret, pseudocount):
    # limits for data input (for testing)
    s=10000000
    e=21000000

    # inpt_ints = inpt_bw.intervals(CHR, s,e)
    for f in inpt_fname:
        inpt_bw=pyBigWig.open(f)
        chrs_in = inpt_bw.chroms().keys()
        assert (set(chrs_in) == set(CHR)), \
            "set of chromosomes in bigwig (%s; %s)\ndiffers from used set (%s)" % \
                    (f, " ".join(chrs_in), " ".join(CHR)) 
        for c in CHR:
            inpt_ints = inpt_bw.intervals(c)
            if ret[c] is None:
                # use pseudocount if required
                init_cnt = np.ones(len(inpt_ints)) if pseudocount else np.zeros(len(inpt_ints))
                ret[c] = {'inpt':init_cnt,
                          'start':np.array([i[0] for i in inpt_ints]),
                          'end':np.array([i[1] for i in inpt_ints]),
                          'header':(c, inpt_bw.chroms(c))}
            ret[c]['inpt']   += np.array([i[2] for i in inpt_ints])
        inpt_bw.close()


def read_cdna(cdna_fname, CHR, ret, lengths):
    # limits for data input (for testing)
    s=10000000
    e=21000000

    # init coverage for cdna
    for c in CHR:
        ret[c]['cdna'] = np.zeros(lengths[c])

    # and read all cdna files and add coverage values
    for fname in cdna_fname:
        cdna_bw=pyBigWig.open(fname)
        # cdna += np.array([i[2] for i in cdna_bw.intervals(CHR, s,e)])[nonzero_idx]
        for c in CHR:
            ret[c]['cdna'] += np.array([i[2] for i in cdna_bw.intervals(c)])
        cdna_bw.close()


def norm_chr(inpt_fname, cdna_fname, CHR, pseudocount=0):
    # function to import and normalize SuRE data (cDNA counts and input counts)

    # import CHR from first input file if CHR is an empty list
    if not CHR:
        inpt_bw=pyBigWig.open(inpt_fname[0])
        chroms = list(inpt_bw.chroms().keys())
        inpt_bw.close()
        for c in chroms:
            CHR.append(c)

    cov = {c:None for c in CHR}

    read_inpt(inpt_fname, CHR, cov, pseudocount)
    read_cdna(cdna_fname, CHR, cov, lengths={c:len(cov[c]['inpt']) for c in CHR})

    # discard data when inpt=0 (if pseudocount is used there are no zeros)
    if pseudocount is 0:
        for c in CHR:
            zero_idx = np.argwhere(cov[c]['inpt']==0)
            cov[c]['inpt']  = np.delete(cov[c]['inpt'], zero_idx)
            cov[c]['start'] = np.delete(cov[c]['start'], zero_idx)
            cov[c]['end']   = np.delete(cov[c]['end'], zero_idx)
            cov[c]['cdna']  = np.delete(cov[c]['cdna'], zero_idx)

    # clean return value, set normalized SuRE value, sum total score
    len_tot = 0.
    cov_tot = 0.
    for c in CHR:
        cov[c]['norm'] = cov[c]['cdna']/cov[c]['inpt']
        del cov[c]['inpt']
        del cov[c]['cdna']
        ll = cov[c]['end'] - cov[c]['start']
        cov_tot += sum(ll * cov[c]['norm'])
        len_tot += sum(ll)

    # compute genome-wide mean score and scale normalized score to mean = 1
    mean = cov_tot/len_tot
    for c in CHR:
        cov[c]['norm'] = cov[c]['norm']/mean

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
                           values=cov[c]["norm"].tolist())
    # close bigwig file
    print("closing bw")
    norm_bw.close()
    print("writing to output finished")


def main(options):
    ## all chromosome names, in correct order
    # CHRS = ['chr'+str(i) for i in (list(range(1,22)) + ['X', 'Y','M'])]
    # CHRS = ['chr'+str(i) for i in (list(range(19,21)))]
    CHR = []

    cov = norm_chr(options.input, options.cdna, CHR, options.pseudo)
    fname = os.path.join(options.out)
    write_bw(cov, fname, CHR)


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))

    options = parse_options()

    main(options)

