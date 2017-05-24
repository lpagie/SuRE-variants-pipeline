# bed2coverage; 'optimized' for SuRE-seq coverage files

## Introduction

In Sure-seq we want to create coverage files such that specified
fragments may have a coverage of 0, but not the entire genome is
necessarilly included in the coverage file. This reflects cases where
certain fragments exist in a SuRE library but are not detected in a
cDNA/plDNA assay.\\
Also, if consecutive fragments have identical coverage they should
nevertheless not be merged into one fragment.\\
As existing coverage methods (R, bedtools, ...) did not implement these
specs I made this tool.

## Usage

### Input

`bed2coverage` takes input from stdin in the form of space/tab separated text.
Each line is expected to contain 4 fields;\\
  * name of chromosome
  * start of fragment
  * end of fragment
  * count (weight) of current fragment

No checking of input is performed; the user is entirely responsible for feeding
correct input.

### Output

`bed2coverage` writes coverage data in bed-like format, using tab separated text. Each line contains 4 fields:\\
  * name of chromosome
  * start of fragment
  * end of fragment
  * coverage

## Method

Computing the coverage is done by initializing a large vector (length 300e6) to
-1. For all input fragments the corresponding elements in the vector are set to
the fragment's weight (if elements were still at initial value) or increased
with the fragment's weight. All unique start/end positions of the fragments are
collected. When all fragments of a chromosome are processed the coverage vector
is split according to start/end positions, and written to stdout.

## Installation

The source file is in `.../usr/local/src/bed2coverage/`. There is also a
Makefile which either makes the executable or installs the executable in
`.../usr/local/bin/`

