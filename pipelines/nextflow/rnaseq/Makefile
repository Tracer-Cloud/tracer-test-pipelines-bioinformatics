

SHELL := bash
SPACK_VERSION := 0.23.0
SPACK_DIR := ../../..

# Get parent of makefile which 
BASEDIR := $(abspath $(CURDIR)/)
CONFIG_DIR := ../config/


.PHONY: download_spack
download_spack:
	@if [ ! -f spack-$(SPACK_VERSION).tar.gz ]; \
	then \
		echo "Downloading Spack $(SPACK_VERSION)"; \
		curl -LOJ https://github.com/spack/spack/releases/download/v$(SPACK_VERSION)/spack-$(SPACK_VERSION).tar.gz; \
	fi; \
	if [ ! -d spack ]; \
	then \
		rm -r spack-$(SPACK_VERSION) || true; \
		tar -xzf spack-$(SPACK_VERSION).tar.gz; \
		mv spack-$(SPACK_VERSION) spack; \
	fi; \
	. spack/share/spack/setup-env.sh; \
	spack bootstrap now

.PHONY: setup_environment
setup_environment: download_spack
	@. $(SPACK_DIR)/share/spack/setup-env.sh; \
	spack env activate -d ..; \
	spack install


.PHONY: test_rnaseq
test_rnaseq:
	@. $(SPACK_DIR)/share/spack/setup-env.sh; \
	spack env activate -d ..; \
	nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		-profile test

.PHONY: test_rnaseq_full
test_rnaseq_full:
	@. $(SPACK_DIR)/share/spack/setup-env.sh; \
	spack env activate -d ..; \
	nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		-profile test_full