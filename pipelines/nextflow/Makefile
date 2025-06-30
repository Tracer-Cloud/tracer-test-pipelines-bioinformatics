SHELL := bash

BASEDIR := $(abspath $(CURDIR)/)
CONFIG_DIR := $(BASEDIR)/config/

.PHONY: check_prerequisites
check_prerequisites:
	@if ! command -v python3 &> /dev/null || ! command -v git &> /dev/null || ! command -v nextflow &> /dev/null; then \
		echo "Error: Python3, Git, and Nextflow are required."; \
		exit 1; \
	fi

.PHONY: test_sarek
test_sarek:
	@nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/sarek \
		-params-file $(CONFIG_DIR)/sarek-params.json \
		--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
		-profile docker,arm,test

.PHONY: test_rnaseq
test_rnaseq:
	@nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
		-profile test

.PHONY: test_rnaseq_full
test_rnaseq_full:
	@nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
		-profile test_full

.PHONY: test_proteinfold
test_proteinfold:
	@nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/proteinfold \
		-params-file $(CONFIG_DIR)/proteinfold-params.json \
		--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
		-profile docker,arm,test

.PHONY: clean
clean:
	@echo "Cleaning up..."
