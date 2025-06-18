# SHELL := bash
SPACK_DIR := $(abspath $(CURDIR)/spack)
ENV_DIR := $(abspath $(CURDIR))
CONFIG_DIR := $(abspath $(CURDIR)/config)

# Get parent of makefile which 
BASEDIR := $(abspath $(CURDIR)/)

.PHONY: check_prerequisites
check_prerequisites:
	@if ! command -v python3 &> /dev/null || ! command -v git &> /dev/null; then \
		echo "Error: Python3 and Git are required."; \
		exit 1; \
	fi

.PHONY: install_spack
install_spack: check_prerequisites
	@if [ ! -d "$(SPACK_DIR)" ]; then \
		echo "Installing Spack..."; \
		git clone --depth=1 https://github.com/spack/spack.git > /dev/null 2>&1; \
		echo "Spack repository cloned successfully."; \
	fi

.PHONY: setup_spack_env
setup_spack_env: install_spack
	@echo "Setting up Spack environment..."
	@if [ ! -f ~/.spack/setup-env.sh ]; then \
		mkdir -p ~/.spack; \
		echo "export SPACK_ROOT=$(SPACK_DIR)" > ~/.spack/setup-env.sh; \
		echo "source \$SPACK_ROOT/share/spack/setup-env.sh" >> ~/.spack/setup-env.sh; \
		echo "export PATH=\$SPACK_ROOT/bin:\$PATH" >> ~/.spack/setup-env.sh; \
		chmod +x ~/.spack/setup-env.sh; \
	fi
	@SHELL_NAME=$(shell basename $$SHELL); \
	if [ "$$SHELL_NAME" = "zsh" ]; then \
		RC_FILE=~/.zshrc; \
	else \
		RC_FILE=~/.bashrc; \
	fi; \
	if ! grep -q "source ~/.spack/setup-env.sh" $$RC_FILE; then \
		echo "source ~/.spack/setup-env.sh" >> $$RC_FILE; \
	fi
	@echo "Verifying Spack installation..."
	@. $(SPACK_DIR)/share/spack/setup-env.sh; \
	if ! command -v spack > /dev/null 2>&1; then \
		echo "Error: 'spack' not found in PATH after sourcing Spack environment."; \
		exit 1; \
	fi; \
	if ! spack --version; then \
		echo "Error: 'spack --version' failed. Spack may not be installed correctly."; \
		exit 1; \
	fi
	@echo "Spack environment setup completed successfully."

.PHONY: verify_spack
verify_spack: setup_spack_env
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && spack --version" > /dev/null 2>&1 || { \
		echo "Error: Spack verification failed"; \
		exit 1; \
	}

.PHONY: create_env
create_env: verify_spack
	@echo "Creating or activating Spack environment..."
	@bash -c 'source "$(SPACK_DIR)/share/spack/setup-env.sh" && \
		cd "$(ENV_DIR)" && \
		if [ -f "spack.yaml" ]; then \
			spack env activate -d . && \
			spack install; \
		else \
			spack env create -d . --with-view "$(ENV_DIR)/.spack-env/view" && \
			spack env activate -d . && \
			spack install; \
		fi'
	@echo "Spack environment ready."

.PHONY: setup_environment
setup_environment: create_env
	@echo "Installing dependencies..."
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		spack install"
	@echo "Environment setup completed successfully."

.PHONY: test_rnaseq
test_rnaseq: setup_environment
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(CONFIG_DIR)/rnaseq-params.json \
			--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
			-profile test"

.PHONY: test_rnaseq_full
test_rnaseq_full: setup_environment
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(CONFIG_DIR)/rnaseq-params.json \
			--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml \
			-profile test_full"

.PHONY: clean
clean:
	@echo "Cleaning up..."
	@rm -rf spack ~/.spack
