SHELL := bash
SPACK_VERSION := 0.23.0

BASEDIR := $(abspath $(CURDIR)/../)
CONFIG_DIR := $(BASEDIR)/config/nextflow

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
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	spack install

.PHONY: test_sarek_aws_batch
test_sarek_aws_batch:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	nextflow -c $(CONFIG_DIR)/batch.config run https://github.com/nf-core/sarek \
		-params-file $(CONFIG_DIR)/sarek-params.json \
		-profile test

.PHONY: test_full_sarek_aws_batch
test_full_sarek_aws_batch:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	nextflow -c $(CONFIG_DIR)/batch.config run https://github.com/nf-core/sarek \
		-params-file $(CONFIG_DIR)/sarek-params.json \
		-profile test_full

.PHONY: test_rnaseq_aws_batch
test_rnaseq_aws_batch:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	nextflow -c $(CONFIG_DIR)/batch.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		-profile test

.PHONY: test_full_rnaseq_aws_batch
test_full_rnaseq_aws_batch:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	nextflow -c $(CONFIG_DIR)/batch.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		-profile test_full

.PHONY: test_proteinfold_aws_batch
test_proteinfold_aws_batch:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	nextflow -c $(CONFIG_DIR)/batch.config -c $(CONFIG_DIR)/proteinfold.config run https://github.com/nf-core/proteinfold \
		-params-file $(CONFIG_DIR)/proteinfold-params.json \
		-profile test \
		--use_gpu

.PHONY: test_full_proteinfold_aws_batch
test_full_proteinfold_aws_batch:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	nextflow -c $(CONFIG_DIR)/batch.config -c $(CONFIG_DIR)/proteinfold.config run https://github.com/nf-core/proteinfold \
		-params-file $(CONFIG_DIR)/proteinfold-params.json \
		-profile test_full \
		--use_gpu
