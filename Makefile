SHELL = /bin/bash
PREFIX = $(HOME)/usr/local

# CONDAROOT = /home/NFS/users/l.pagie/usr/local/src/miniconda2
CONDAROOT = /DATA/home/ludo/miniconda3/

.PHONY: install
install: subdirs SuRE-snakemake
	sed -i.org 's~^\(CODE_BASE\s\+=\s\+\).*~\1"'"$${PWD}/code"'"~g' SuRE-snakemake
	source $(CONDAROOT)/bin/activate && conda env create -n SuRE-pipeline -f code/conda-wasp-environment.yml && source deactivate

SUBDIRS = code/bed2coverage

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS): 
	$(MAKE) -C $@

