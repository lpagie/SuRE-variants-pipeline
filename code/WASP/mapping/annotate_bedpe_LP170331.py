# AUTHOR/DATE

# INTRO/BACKGROUND
#   A script to annotate a SuRE-seq combined bepe file with SNP variant calls.
#   The input bedpe file contains the sequences of the forw/rev reads which
#   define the SuRE fragment. For the entire fragment SNP positions are
#   determined and the known allelic variants are compared to the read
#   sequences. The SNP positions (relative to start of the SuRE frasgment,
#   0-based) and observed base identities, are added to the input bedpe file.
#   Also added are flags for each SNP position whether that position is in the
#   (forw/rev/both) reads or not, and a number indicating the concordance of
#   base identities to an inferred maternal/paternal chromosome.

# USAGE / INPUT / ARGUMENTS / OUTPUT 
# USAGE:
# INPUT:
#  iPCR-combined-bedpe file: tab separated text file. each row contains:
#    1 chr
#    2 start
#    3 end
#    4 length
#    5 strand
#    6 barcode
#    7 count
#    8 internal end
#    9 internal start
#    10 MAPQ
#    11 MD1
#    12 MD2
#    13 T/F (is there an alternative alignment possible?)
#    14 T/F (is there an alternative alignment possible?)
#    15 seq1 (always plus strand, reverse direction)
#    16 seq2 (always minus strand, reverse direction?)
# ARGUMENTS:
# OUTPUT:

# VERSIONS:
# TODO:
#   - adapt def write_results (but probably results are written during the read loop, in which case this function is superfluous))
#   - 
#   - 
#   - 
#   - 

import sys
import argparse
import numpy as np
# import pysam

import util
import snptable

import tables
import gzip

import os
import re
import string

def write_results(out_f, chrom_name, snp_tab, ref_matches,
                  alt_matches, oth_matches, geno_sample):

    haps = None
    has_haps = False
    
    if geno_sample:
        # get index for this sample in the haplotype table
        samp_idx_dict = dict(zip(snp_tab.samples,
                                 range(len(snp_tab.samples))))

        if geno_sample in samp_idx_dict:
            idx = samp_idx_dict[geno_sample]
            geno_hap_idx = np.array([idx*2, idx*2+1], dtype=np.int)
            haps = snp_tab.haplotypes[:,geno_hap_idx]
            has_haps = True
            sys.stderr.write("geno_hap_idx: %s\n" % repr(geno_hap_idx))
        else:
            sys.stderr.write("WARNING: sample %s is not present for "
                             "chromosome %s\n" % (geno_sample, chrom_name))
            haps = None
            has_haps = False

    for i in range(snp_tab.n_snp):
        if has_haps:
            geno_str = "%d|%d" % (haps[i, 0], haps[i, 1])
        else:
            geno_str = "NA"
        out_f.write("%s %d %s %s %s %d %d %d\n" %
                    (chrom_name, snp_tab.snp_pos[i],
                     snp_tab.snp_allele1[i], snp_tab.snp_allele2[i],
                     geno_str, ref_matches[i], alt_matches[i],
                     oth_matches[i]))


def write_header(out_f):
    out_f.write("CHROM SNP.POS REF.ALLELE ALT.ALLELE GENOTYPE REF.COUNT "
                "ALT.COUNT OTHER.COUNT\n")


    
def parse_samples(samples_str):
    """Gets list of samples from --samples argument. This may be 
    a comma-delimited string or a path to a file. If a file is provided 
    then the first column of the file is assumed to be the sample name"""

    if samples_str is None:
        return None
        
    # first check if this is a path to a file
    if os.path.exists(samples_str) and not os.path.isdir(samples_str):
        samples = []

        if util.is_gzipped(samples_str):
            f = gzip.open(samples_str)
        else:
            f = open(samples_str)

        for line in f:
            # assume first token in line is sample name
            samples.append(line.split()[0])

        sys.stderr.write("read %d sample names from file '%s'\n" %
                         (len(samples), samples_str))
                    
        f.close()
    else:    
        # otherwise assume comma-delimited string
        if ("," not in samples_str and len(samples_str) > 15) \
           or ("/" in samples_str):
            sys.stderr.write("WARNING: --samples argument (%s) "
                             "does not look like sample name "
                             "but is not path to valid file. "
                             "Assuming it is a sample name anyway."
                             % samples_str)

        samples = samples_str.split(",")
        sys.stderr.write("SAMPLES: %s\n"% repr(samples))


    return samples


    


def parse_options():
    parser = argparse.ArgumentParser(description=" This script takes a bedpe-like text filei from SuRE-seq as input"
        "and annotates the SuRE fragments with SNP variants present in the "
        "fragment.  Output is written to stdout.")

    parser.add_argument("--snp_dir", action='store', 
                        help=("Directory containing SNP text files "
                              "This directory should contain one file per "
                              "chromosome named like chr<#>.snps.txt.gz. "
                              "Each file should contain 3 columns: position "
                              "RefAllele AltAllele"),
                        default=None)
        

    parser.add_argument("--snp_tab",
                        help="Path to HDF5 file to read SNP information "
                        "from. Each row of SNP table contains SNP name "
                        "(rs_id), position, allele1, allele2.",
                        metavar="SNP_TABLE_H5_FILE",
                        default=None)
    
    parser.add_argument("--snp_index",
                        help="Path to HDF5 file containing SNP index. The "
                        "SNP index is used to convert the genomic position "
                        "of a SNP to its corresponding row in the haplotype "
                        "and snp_tab HDF5 files.",
                        metavar="SNP_INDEX_H5_FILE",
                        default=None)
    
    parser.add_argument("--haplotype",
                        help="Path to HDF5 file to read phased haplotypes "
                        "from. When generating alternative reads "
                        "use known haplotypes from this file rather "
                        "than all possible allelic combinations.",
                        metavar="HAPLOTYPE_H5_FILE",
                        default=None)

    parser.add_argument("--samples",
                        help="Use only haplotypes and SNPs that are "
                        "polymorphic in these samples. "
                        "SAMPLES can either be a comma-delimited string "
                        "of sample names or a path to a file with one sample "
                        "name per line (file is assumed to be "
                        "whitespace-delimited and first column is assumed to "
                        "be sample name). Sample names should match those "
                        "present in the haplotype HDF5 file. Samples are "
                        "ignored if no haplotype file is provided.",
                        metavar="SAMPLES", default=None)


    parser.add_argument("--genotype_sample",
                        metavar="GENO_SAMPLE",
                        help="output genotypes for sample with name "
                        "GENO_SAMPLE alongside allele-specific counts. "
                        "GENO_SAMPLE must match one "
                        "of the names present in the haplotype HDF5 file. "
                        "If the --samples argument is provided then "
                        "GENO_SAMPLE must also be one of the specified "
                        "samples. If --genotype_sample is "
                        "not provided or the GENO_SAMPLE does not match any "
                        "of the samples in haplotype file then NA is "
                        "output for genotype.", default=None)
        
    parser.add_argument("bam_filename", action='store',
                        help="Coordinate-sorted input BAM file "
                        "containing mapped reads.")

#    parser.add_argument("bedpe-combined_filename", action='store',
#                        help="Coordinate-sorted input bedpe-combined file "
#                        "containing SuRE fragments.")


    options = parser.parse_args()
    
    if options.snp_dir:
        if(options.snp_tab or options.snp_index or options.haplotype):
            parser.error("expected --snp_dir OR (--snp_tab, --snp_index and "
                         "--haplotype) arguments but not both")
    else:
        if not (options.snp_tab and options.snp_index and options.haplotype):
            parser.error("either --snp_dir OR (--snp_tab, "
                         "--snp_index AND --haplotype) arguments must be "
                         "provided")
     
    return options


def cigar2tuple(cigar):
    # returned cigar is a list of tuples. Each tuple has two entries. The first
    # entry specifies the character in the cigar and the second entry
    # specifies the length of that character. The values are
    # M       BAM_CMATCH      0
    # I       BAM_CINS        1
    # D       BAM_CDEL        2
    # N       BAM_CREF_SKIP   3
    # S       BAM_CSOFT_CLIP  4
    # H       BAM_CHARD_CLIP  5
    # P       BAM_CPAD        6
    # =       BAM_CEQUAL      7
    # X       BAM_CDIFF       8
    # E.g. (0, 5) means 5 matches, and (4, 2) means a soft clip of 2bp
    #
    # Map from CIGAR symbols to corresponding BAM values
    symbolmap = {s:l for s,l in zip(('M','I','D','N','S','H','P','=','X'), (0,1,2,3,4,5,6,7,8))}
    #
    # regexp for splitting cigar string in substrings "LENGTH|OPERATOR"
    split_re='\d+[MIDNSHP=X]'
    #
    # iterate over list of subcigar strings, coerce lengths into int and operator into BAM value
    return [(symbolmap[i[-1]],int(i[0:-1])) for i in re.findall(split_re, cigar)]

def iupac(base1, base2):
    if (base1 == base2): return base1
    if ((base1=="A" and base2=="G") or (base2=="A" and base1=="G")): return "R"
    if ((base1=="A" and base2=="T") or (base2=="A" and base1=="T")): return "W"
    if ((base1=="A" and base2=="C") or (base2=="A" and base1=="C")): return "M"
    if ((base1=="G" and base2=="C") or (base2=="G" and base1=="C")): return "S"
    if ((base1=="G" and base2=="T") or (base2=="G" and base1=="T")): return "K"
    if ((base1=="C" and base2=="T") or (base2=="C" and base1=="T")): return "K"

trans = string.maketrans('ATGC', 'TACG')
def reverse_complement(dna):
    # complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # return ''.join([complement[base] for base in dna[::-1]])
    return(dna[::-1].translate(trans))

def print_comma_sep_list(lst):
    #
    if(len(lst) == 0):
        return ""
    #
    string=str(lst[0])
    for e in lst[1:]:
        string +=","+str(e)
    #
    return(string)

def main(bedpe_filename, snp_dir=None):

    endl = os.linesep

    out_f = sys.stdout
    
    # bam = pysam.Samfile(bam_filename)
    bedpe = gzip.open(bedpe_filename)
        
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])
    snp_chrom = None

    snp_tab = snptable.SNPTable()
    read_pair_cache = {}

    for line in bedpe:
        # print line

        line = line.rstrip(endl)
        cur_chrom,c_start,c_end,_,c_strand,_,_,c_iend,c_istart,_,_,_,_,_,c_seq1,c_seq2,c_cigar1,c_cigar2 = line.split("\t")
        # cur_chrom,c_start,c_end,_,c_strand,_,_,c_iend,c_istart,_,_,_,_,_,c_seq1,c_seq2 = line.split("\t")

        c_start=int(c_start)
        c_istart=int(c_istart)
        c_end=int(c_end)
        c_iend=int(c_iend)

        if (len(c_seq1) != (c_iend-c_start+1) or len(c_seq2) != (c_end-c_istart+1)):
            # print "indels in "+line
            continue
        c_cigar1=str(len(c_seq1))+"M"
        c_cigar2=str(len(c_seq2))+"M"

        # one of the reads is on minus strand, depending on the strand to which
        # the fragment is mapped; I need all reads on the plus strand
        c_seq2=reverse_complement(c_seq2)
        # print("SEQ = "+c_seq2)

        if (snp_chrom is None) or (cur_chrom != snp_chrom):
            # this is a new chromosome
            
            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            # cur_tid = read.tid
            snp_chrom = cur_chrom
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)

            # read SNPs for next chromomsome
            # read SNPs from text file
            snp_filename = "%s/%s.snps.txt.gz" % (snp_dir, cur_chrom)
            snp_tab.read_file(snp_filename)

            sys.stderr.write("read %d SNPs\n" % snp_tab.n_snp)
            
        # loop over all SNP that overlap this read; record:
        # - read_pos: SNP position in current read (ie read1 or read2), used to determine base identities
        # - frag_pos: SNP position in entire SuRE-fragment, may occur multiple times if both reads overlap same SNP
        # - read_base: base identities in current read, is compared to allele-variants for SNPs in current read
        # - frag_base: base identities for all SNPs in SuRE-fragment
        # - snp_pos: chromosome position of SNPs
        # - snp_var: whether frag_base is reference allele (0), alternative allele (1), or non-matching (2), or unknown (3)
        snp_idx, snp_read_pos, indel_idx, indel_read_pos = \
                snp_tab.get_overlapping_snps_from_bedpe(c_start-1, cigar2tuple(c_cigar1), len(c_seq1))
        read_pos = [p-1 for p in snp_read_pos]
        frag_pos = read_pos
        read_base = [c_seq1[p] for p in read_pos]
        frag_base = read_base
        snp_pos = [snp_tab.snp_pos[i] for i in snp_idx]
        snp_var = [int((b==snp_tab.snp_allele1[i] and '0') or (b==snp_tab.snp_allele2[i] and 1) or 2) for b, i in zip(read_base, snp_idx)]
        snp_ind = snp_idx

        snp_idx, snp_read_pos, indel_idx, indel_read_pos = \
                snp_tab.get_overlapping_snps_from_bedpe(c_istart-1, cigar2tuple(c_cigar2), len(c_seq2))
        read_pos = [p-1 for p in snp_read_pos]
        frag_pos = frag_pos + [p+(c_istart - c_start) for p in read_pos]
        read_base = [c_seq2[p] for p in read_pos]
        frag_base = frag_base + read_base
        snp_pos = snp_pos + [snp_tab.snp_pos[i] for i in snp_idx]
        snp_var = snp_var + \
                [int((b==snp_tab.snp_allele1[i] and '0') or (b==snp_tab.snp_allele2[i] and 1) or 2) for b, i in zip(read_base, snp_idx)]
        snp_ind += snp_idx

        # if reads do not overlap the sequence in between is also checked for SNP positions
        if c_iend < (c_istart-1):
            l = c_istart - c_iend - 1
            c = str(l)+"M"
            snp_idx, snp_read_pos, indel_idx, indel_read_pos = \
                    snp_tab.get_overlapping_snps_from_bedpe(c_iend+1, cigar2tuple(c), l)
            frag_pos = frag_pos + [p+c_iend-c_start+1 for p in snp_read_pos]
            frag_base = frag_base + [iupac(snp_tab.snp_allele1[i], snp_tab.snp_allele2[i]) for i in snp_idx]
            # frag_base = frag_base + [snp_tab.snp_allele1[i]+snp_tab.snp_allele2[i] for i in snp_idx]
            snp_pos = snp_pos + [snp_tab.snp_pos[i] for i in snp_idx]
            snp_var = snp_var + [3 for  i in snp_idx]
            snp_ind += snp_idx

        line = line + "\t"+ print_comma_sep_list(frag_pos)+"\t"+ print_comma_sep_list(frag_base)+"\t"+ print_comma_sep_list(snp_pos)+"\t"+ print_comma_sep_list(snp_var)+"\t"+print_comma_sep_list(snp_ind)
        # line = line + "\t"+ print_comma_sep_list(frag_pos)+"\t"+ print_comma_sep_list(frag_base)+"\t"+ print_comma_sep_list(snp_var)
        print(line)
#        print(line)
#        # frag_pos
#        print("\t"+ print_comma_sep_list(frag_pos))
#        # frag_base
#        print("\t"+ print_comma_sep_list(frag_base))
#        # snp_pos
#        print("\t"+ print_comma_sep_list(snp_pos))
#        # snp_var
#        print("\t"+ print_comma_sep_list(snp_var))
#        print("\n")


if __name__ == "__main__":
    sys.stderr.write("command: %s\n" % " ".join(sys.argv))

    options = parse_options()
    samples = parse_samples(options.samples)

    
#    main(options.bam_filename, 
#         snp_dir=options.snp_dir,
#         snp_tab_filename=options.snp_tab,
#         snp_index_filename=options.snp_index,
#         haplotype_filename=options.haplotype,
#         samples=samples, geno_sample=options.genotype_sample)
    
    main(options.bam_filename, snp_dir=options.snp_dir)
    

    
