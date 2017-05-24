/* Ludo Pagie, 170116, bed2coverage.c
 *
 * FUNCTIONALITY:
 *  A program which computes a coverage file (bedgraph format) from a set of
 *  (overlapping) genomic fragments with associated counts (possibly zero) in
 *  bed format.
 *
 *  The use case is input from SuRE data. I.e. the input fragments are iPCR
 *  fragments. The resulting output regions are the minimal regions defined by a
 *  unique set of SuRE fragments overlapping with this region. The input
 *  fragments are named 'fragment', and the output regions are named 'region'. 
 *
 *  The input fragments are expected to be 1-based, fully-closed. The output
 *  regions are 0-based and half-open
 *  (http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/).
 *  This means that of two consecutive regions the end position of the first
 *  coincides with the start position of the second. Transitions from one region
 *  to the next are defined as start- and end- points of SuRE-fragments. Such
 *  transition points may coincide (possibly multiple times). The start of a
 *  fragment at genomic position x (1-based coordinate) defines the start of a
 *  region at position x-1 (0-based coordinate). The same position x-1 defines
 *  the end of the preceding region, as the coordinate system is half-open.
 *
 * USAGE
 *  cat inputfile | bed2wig > outputfile
 *
 * INPUT
 *  Input is read from stdin. The input should be space or tab separated text, with 4 fields per line:
 *   - name of space (chromosome)
 *   - start position
 *   - end position
 *   - count
 *  genomic positions should be using a 1-based fully-closed coordinate system.
 *  The program does NOT allow header lines
 *
 * OUTPUT
 *  Output is written on stdout, in bedGraph format. I.e.
 *  chromosome/start/end/count, tab separated.
 *
 *  NOT ANYMORE!!
 *  ! Each data section corresponding to a chromosome (space) starts with a (kind
 *  ! of?) definition line (eg "#bedGraph section chr1:61911-249237851"). This
 *  ! format should be valid input for the UCSC tool wigToBigWig
 *  NOT ANYMORE!!
 *
 * METHOD
 *  Simply collect counts for all SuRE fragments in a single vector and split
 *  it into regions according to start/end positions of input fragments.
 *
 * VERSIONS
 *  0.1
 *
 * TODO
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

const int maxlen=300000000; // max length of an input chromosome, ie the 'space' on which fragments are placed

int intcmp(const void *aa, const void *bb)
{
      const int *a = aa, *b = bb;
          return (*a < *b) ? -1 : (*a > *b);
}

void writeRegionsInSpace(int *regionAll, int regionCount, int *regionUniq, int *cov, char *lab) {
/* Function to write all regions defined by their respective end positions in
 * the array 'regionAll' (after removing duplicates), and their coverage (array
 * 'cov') >= 0
 *
 */
  int i, j;
  int start, end;
  int regionInd=0;

  // sort the array with region end positions
  qsort(regionAll, regionCount, sizeof(int), intcmp);
  // remove duplicates; copy the uniq elements to array regionUniq
  regionUniq[0]=regionAll[0];
  for (i=1; i<regionCount; i++)
    if (regionAll[i]>regionUniq[regionInd])
      regionUniq[++regionInd]=regionAll[i];
  // fprintf(stderr, "after uniq %d items left\n", regionInd);

  // print bed header
  // printf("#bedGraph section %s:%d-%d\n", lab, regionUniq[0], regionUniq[regionInd-1]);

  // iterate over all end positions and print the region (if coverage > -1)
  // As the regions are 0-based and half-open consucutive regions end and start at the exact same position
  start=regionUniq[0];
  for(i=1; i<regionInd; i++) {
    end=regionUniq[i];
    if (cov[end] > -1) {
      printf("%s\t%d\t%d\t%d\n", lab, start, end, cov[end]);
    }
    // reset coverage for this region to prevent this region to be printed
    // multiple times (and to initialize memory for next chromosome)
    for (j=start; j<=end;j++)
      cov[j]=-1;
    start=end;
  }
  // fprintf(stderr, "for chr %s wrote %d fragments\n", lab, i);
}
 
int main() 
{
  int i, ret, regionInd=0, ind;                 // local var's, counter, etc
  int *cov=malloc(maxlen*sizeof(int));          // vector with coverage count
  int *regionAll=malloc(maxlen*sizeof(int));    // vector for collecting start/end positions of all regions
  int *regionbuffer=malloc(maxlen*sizeof(int)); // memory buffer for function writeRegionsInSpace
  int start, end, cnt, cc=0;                    // variables for reading input data
  char CHR[1000], CHRprev[1000];                // variables for reading input data

  // initialize coverage vector to '-1', indicating positions not covered by any fragment
  for(i=0; i<maxlen; i++)
      cov[i]=-1;
  // initialize CHR/CHRprev to enable testing for reading new 'space label'
  memset(CHR, '\0', sizeof(CHR));
  memset(CHRprev, '\0', sizeof(CHRprev));

  // read input, update coverage vector, record transition points of (overlap of) SuRE fragments
  while((ret=scanf("%s%d%d%d\n",CHR, &start, &end, &cnt))!=EOF) {
    cc++;
    // bad input
    if(ret!=4) {
      fprintf(stderr, "foutje\n");
      fprintf(stderr, "%s %d %d %d %d\n", CHR, start, end, cnt, ret);
      return -1;
    }
    if (strcmp(CHR, CHRprev) !=0) {
      //fprintf(stderr, "new chr=%s, old chr=%s, at input line %d\n", CHR, CHRprev, cc);
      if (strlen(CHRprev)>0) {
        // check whether we see a new 'space label'. If so, write the data for the previous 'space' and reset iterators etc
        writeRegionsInSpace(regionAll, regionInd, regionbuffer, cov, CHRprev);
        regionInd=0;
	// initialize coverage vector to '-1', indicating positions not covered by any fragment
	for(i=0; i<maxlen; i++)
	  cov[i]=-1;
      }
      strcpy(CHRprev, CHR);
    }
    // update coverage vector
    for(i=start; i<=end; i++) {
      if(cov[i] == -1) 
        cov[i]=cnt;
      else
        cov[i] += cnt;
    }
    // record the start/end positions of this fragment (Nb, the start-1 is
    // needed to convert to 0-based coordinates, the end position remains the
    // same as the coordinate system for regions is half-open)
    regionAll[regionInd++]=start-1;
    regionAll[regionInd++]=end;
  }
  writeRegionsInSpace(regionAll, regionInd, regionbuffer, cov, CHRprev);

  return 0;
}
