.PHONY: install
install: subdirs SuRE-snakemake

# sed -i.org 's~^\(CODE_BASE\s\+=\s\+\).*~\1"'"$${PWD}/code"'"~g' SuRE-snakemake

SUBDIRS = code/bed2coverage

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS): 
	$(MAKE) -C $@

check:
	@type conda >/dev/null 2>&1|| { \
		echo "Conda is required to build this application"; \
		exit 1; \
	}
	conda env update -n SuRE-pipeline -f code/conda-env-sure.yml 
	conda env update -n SuRE-pipeline -f code/conda-env-py2.yml 
