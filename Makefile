PREFIX = $(HOME)/usr/local

.PHONY: install
install: subdirs SuRE-snakemake
	sed -i.org 's~^\(CODE_BASE\s\+=\s\+\).*~\1"'"$${PWD}/code"'"~g' SuRE-snakemake

SUBDIRS = code/bed2coverage

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS): 
	$(MAKE) -C $@

