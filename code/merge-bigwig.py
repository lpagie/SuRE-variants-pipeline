import sys
import pyBigWig
import argparse
import os.path
import natsort

def parse_options():
    # parse user options:

    parser = argparse.ArgumentParser(
        description="merge a set of bigwig files containing data for different chromosomes")
    parser.add_argument("--inbw", 
                        required=True,
                        nargs='+', 
                        help=("(multiple) filenames of bigwig files for merging"))

    parser.add_argument("--outfname", 
                        required=True,
                        help=("Name for output bigwig file")) 

    options = parser.parse_args()
    return options


def main(options):
    # open all files
    # read all headers
    fnames = options.inbw
    bw     = {}
    header = {}
    chroms = {}
    for fname in fnames:
        print("opening file for input: "+os.path.split(fname)[1])
        bw[fname]     = pyBigWig.open(fname)
        header[fname] = bw[fname].chroms()
        chroms[fname] = list(header[fname].keys())[0]

    # define order based on chromosome names extracted from headers
    idx = natsort.index_natsorted(chroms)
    fnames = natsort.order_by_index(fnames, idx)

    # open bigwig for output
    print("opening bigwig file for output (%s)" % options.outfname)
    out_bw = pyBigWig.open(options.outfname, "w")
    assert(out_bw is not None)
    # construct sorted header and write to output
    header = [list(header[f].items())[0] for f in fnames]
    print(str(header))
    out_bw.addHeader(header)

    # loop over sorted chromosome-names/file-names
    ## import data from input bw
    ## add read data to output bw
    for fname in fnames:
        print("exporting data from chrom; "+chroms[fname]+" (file: "+fname+")")
        ints = bw[fname].intervals(chroms[fname])
        chrs = [chroms[fname]]*len(ints)
        out_bw.addEntries(chrs, [i[0] for i in ints], ends=[i[1] for i in ints], values=[i[2] for i in ints])
        # st, en, sc = zip(*ints)
        # out_bw.addEntries(chrs, list(st), list(en), list(sc))
        # out_bw.addEntries(chrs, *zip(*ints))

    print("closing bigwig file for output (%s)" % options.outfname)
    out_bw.close()
    return True

if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))

    options = parse_options()

    main(options)


