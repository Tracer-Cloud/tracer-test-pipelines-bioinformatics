# SHELL := bash
ENV_DIR := $(abspath $(CURDIR))
CONFIG_DIR := $(abspath $(CURDIR)/config)
CONDA_ENV_NAME := nextflow-minimal
CONDA_ENV_FILE := $(abspath $(CURDIR)/environment.yml)
CONDA_BIN := $(HOME)/miniconda3/bin/conda

.PHONY: check_prerequisites
check_prerequisites:
	@if ! command -v python3 &> /dev/null || ! command -v git &> /dev/null; then \
		echo "Error: Python3 and Git are required."; \
		exit 1; \
	fi

.PHONY: install_conda
install_conda:
	@if [ ! -x "$(CONDA_BIN)" ]; then \
		echo "Conda not found. Installing Miniconda..."; \
		bash ../macos-arm64/nextflow-conda/setup.sh; \
	else \
		echo "Conda already installed at $(CONDA_BIN)"; \
	fi

.PHONY: create_conda_env
create_conda_env: install_conda
	@echo "Creating or updating conda environment..."
	@$(CONDA_BIN) env update -f $(CONDA_ENV_FILE) --name $(CONDA_ENV_NAME) --prune
	@echo "Conda environment '$(CONDA_ENV_NAME)' is ready."

.PHONY: setup_environment
setup_environment: create_conda_env
	@echo "Environment setup completed successfully."

.PHONY: test_rnaseq
test_rnaseq: setup_environment
	@bash -c "source $(HOME)/miniconda3/etc/profile.d/conda.sh && conda activate $(CONDA_ENV_NAME) && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
		-profile test"

.PHONY: test_rnaseq_full
test_rnaseq_full: setup_environment
	@bash -c "source $(HOME)/miniconda3/etc/profile.d/conda.sh && conda activate $(CONDA_ENV_NAME) && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
		-params-file $(CONFIG_DIR)/rnaseq-params.json \
		--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
		-profile test_full"

.PHONY: clean
clean:
	@echo "Cleaning up..."
	@$(CONDA_BIN) env remove -n $(CONDA_ENV_NAME) || true
