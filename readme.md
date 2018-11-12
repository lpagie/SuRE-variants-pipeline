<!--pandoc
t: html
toc:
s:
self-contained:
highlight-style:tango
-->

# (former) SUBMODULE ALERT!!!!!

This pipeline uses WASP (https://github.com/bmvdgeijn/WASP) for some of the SNP
related data processing. I have a forked/adapted version of WASP
(https://github.com/lpagie/WASP) which I added to this repos first as a submodule.

As the submodule approach often confused me I opted to simplify it by including
the code directly into this repos. I used the zip archive of my cloned
WASP-repos
(https://github.com/lpagie/WASP/commit/fe5113d5077d66d055c9a97943e08fcd7881ab24),
unpacked it into `code/WASP`, and added it to this repos. i


# DISCLAIMER ALERT!!!

This collection of scripts is far from a production level set of software.
Developed over a few years it got the job done. Many times along the way I lost
my focus wrt proper documentation and other proper-coding methods. Also, as a
pipeline build over the years with increasing required functionalities, etc, it
is far from a computational optimal pipeline as well. For instance, I wrote my
first python code in this pipeline.

As a result it is hardly to be expected that this pipeline will directly
function as intended after cloning the repos. I will do my best though to
assist in case you want to use it for your work.


# Download and install

In the process of developing this pipeline I included the WASP repos. This
included quite a large history, resulting in a large repos (>300Mb). The
actuale source tree is a bit over 1Mb. As far as I know there are two ways to
limit the doenload size; 1) download only the source code in a zip archive
(<https://github.com/lpagie/SuRE-pipeline/archive/snakemake.zip>), or 2) clone
a *shallow* version of the repos to also get the commit history, etc: `git
clone --depth 4 repos-URL`.

## Compile

After downloading or cloning the repos you need to run `make` in order to
compile one binary in the pipeline:
```
cd toplevel-repos-dir
make install
```

# Overview of processing steps

For an overview of the pipeline see the wiki
