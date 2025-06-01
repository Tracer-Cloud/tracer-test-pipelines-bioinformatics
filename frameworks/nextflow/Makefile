SHELL := bash
SPACK_VERSION := 0.23.0

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

.PHONY: test_sarek
test_sarek:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	cd frameworks/nextflow && nextflow -c ../config/nextflow/local.config run ../pipelines/nf-core/sarek/main.nf \
		-params-file ../config/nextflow/sarek-params.json \
		-profile docker,arm,test


.PHONY: test_rnaseq
test_rnaseq:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	cd frameworks/nextflow && nextflow -c ../config/nextflow/local.config run ../pipelines/nf-core/rnaseq/main.nf \
		-params-file ../config/nextflow/rnaseq-params.json \
		-profile test

.PHONY: test_rnaseq_full
test_rnaseq_full:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	cd frameworks/nextflow && nextflow -c ../config/nextflow/local.config run ../pipelines/nf-core/rnaseq/main.nf \
		-params-file ../config/nextflow/rnaseq-params.json \
		-profile test_full


.PHONY: test_proteinfold
test_proteinfold:
	@. spack/share/spack/setup-env.sh; \
	spack env activate -d .; \
	cd frameworks/nextflow && nextflow -c ../config/nextflow/local.config run ../pipelines/nf-core/proteinfold/main.nf \
		-params-file ../config/nextflow/proteinfold-params.json \
		-profile docker,arm,test
