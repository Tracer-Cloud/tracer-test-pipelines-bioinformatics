.PHONY: all test_rnaseq_light test_rnaseq_quant_only test_rnaseq_align_count setup_environment check_deps spack_env

all: test_rnaseq_light test_rnaseq_quant_only test_rnaseq_align_count

REPO_ROOT := $(abspath ../..)
SPACK_ROOT := $(REPO_ROOT)/spack

check_deps:
	@if ! command -v nextflow > /dev/null 2>&1 || ! command -v java > /dev/null 2>&1 || ! command -v spack > /dev/null 2>&1; then \
		echo "Installing missing dependencies with dependencies_installation.sh..."; \
		bash ../../dependencies_installation.sh; \
	fi

spack_env: check_deps
	@echo "Activating and installing Spack environment..."
	@export SPACK_ROOT=$(SPACK_ROOT); \
	source $$SPACK_ROOT/share/spack/setup-env.sh; \
	if [ -f "spack.yaml" ]; then \
		spack env activate -d . && spack install; \
	else \
		spack env create -d . --with-view "./.spack-env/view" && spack env activate -d . && spack install; \
	fi
	@echo "Spack environment ready."

setup_environment: spack_env

# Run lightweight RNA-seq test pipeline
.PHONY: test_rnaseq_light
test_rnaseq_light: setup_environment
	@export SPACK_ROOT=$(SPACK_ROOT); \
	source $$SPACK_ROOT/share/spack/setup-env.sh; \
	spack env activate -d .; \
	cd test-rnaseq-light && nextflow run main.nf

# Run quantification-only RNA-seq pipeline
.PHONY: test_rnaseq_quant_only
test_rnaseq_quant_only: setup_environment
	@export SPACK_ROOT=$(SPACK_ROOT); \
	source $$SPACK_ROOT/share/spack/setup-env.sh; \
	spack env activate -d .; \
	cd test-rnaseq-quant-only && nextflow run main.nf

# Run alignment + quantification RNA-seq pipeline
.PHONY: test_rnaseq_align_count
test_rnaseq_align_count: setup_environment
	@export SPACK_ROOT=$(SPACK_ROOT); \
	source $$SPACK_ROOT/share/spack/setup-env.sh; \
	spack env activate -d .; \
	cd test-rnaseq-align-count && nextflow run main.nf