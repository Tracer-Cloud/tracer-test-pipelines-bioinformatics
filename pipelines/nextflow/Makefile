SHELL := bash
ENV_DIR := $(abspath $(CURDIR)/..)
CONFIG_DIR := $(abspath $(ENV_DIR)/config)

.PHONY: check_prerequisites
check_prerequisites:
	@if ! command -v python3 &> /dev/null || ! command -v git &> /dev/null || ! command -v nextflow &> /dev/null; then \
		echo "Error: Python3, Git, and Nextflow are required."; \
		exit 1; \
	fi

.PHONY: test_rnaseq
test_rnaseq: check_prerequisites
	@nextflow run https://github.com/nf-core/rnaseq \
		-profile test \
		--outdir results

.PHONY: test_rnaseq_full
test_rnaseq_full: check_prerequisites
	@nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		-profile test_full

.PHONY: clean
clean:
	@echo "Cleaning up..."
