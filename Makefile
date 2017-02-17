PREFIX = $(HOME)/usr/local

install: SuRE-snakemake
	sed -i.org 's~^\(CODE_BASE\s\+=\s\+\).*~\1"'"$${PWD}/code"'"~g' SuRE-snakemake
