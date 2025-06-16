SHELL := bash
SPACK_DIR := $(abspath $(CURDIR)/spack)
ROOT_DIR := $(abspath $(CURDIR)/../../..)
ENV_DIR := $(ROOT_DIR)
CONFIG_DIR := $(abspath $(CURDIR)/../config)

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
		echo "source \$$SPACK_ROOT/share/spack/setup-env.sh" >> ~/.spack/setup-env.sh; \
		echo "export PATH=\$$SPACK_ROOT/bin:\$$PATH" >> ~/.spack/setup-env.sh; \
		chmod +x ~/.spack/setup-env.sh; \
	fi
	@if [ -f ~/.bashrc ]; then \
		if ! grep -q "source.*spack/setup-env.sh" ~/.bashrc; then \
			echo "source ~/.spack/setup-env.sh" >> ~/.bashrc; \
		fi \
	fi
	@if [ -f ~/.zshrc ]; then \
		if ! grep -q "source.*spack/setup-env.sh" ~/.zshrc; then \
			echo "source ~/.spack/setup-env.sh" >> ~/.zshrc; \
		fi \
	fi
	@echo "Verifying Spack installation..."
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && spack --version" > /dev/null 2>&1 || { \
		echo "Error: Failed to set up Spack environment"; \
		exit 1; \
	}
	@echo "Spack environment setup completed successfully."

.PHONY: verify_spack
verify_spack: setup_spack_env
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && spack --version" > /dev/null 2>&1 || { \
		echo "Error: Spack verification failed"; \
		exit 1; \
	}

.PHONY: create_env
create_env: verify_spack
	@echo "Creating Spack environment..."
	@bash -c 'source "$(SPACK_DIR)/share/spack/setup-env.sh" && \
		cd "$(ENV_DIR)" && \
		if [ -d ".spack-env" ]; then \
			rm -rf .spack-env; \
		fi && \
		if [ -f "spack.yaml" ]; then \
			mv spack.yaml spack.yaml.bak; \
		fi && \
		if [ -f "spack.lock" ]; then \
			mv spack.lock spack.lock.bak; \
		fi && \
		spack env create -d . && \
		spack env activate -d . && \
		spack add nextflow@24.10.0 python@3.12 py-pip && \
		spack concretize -f && \
		spack install; \
		if [ -f "spack.yaml.bak" ]; then \
			mv spack.yaml.bak spack.yaml; \
		fi && \
		if [ -f "spack.lock.bak" ]; then \
			mv spack.lock.bak spack.lock; \
		fi'
	@echo "Spack environment created successfully."

.PHONY: setup_environment
setup_environment: create_env
	@echo "Installing dependencies..."
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		spack concretize -f && \
		spack install"
	@echo "Environment setup completed successfully."

.PHONY: test_rnaseq
test_rnaseq: setup_environment
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(BASEDIR) && \
		spack env activate -d $(ENV_DIR) && \
		chmod -R u+w $(BASEDIR)/work && \
		nextflow -c $(CONFIG_DIR)/local.config -c $(CONFIG_DIR)/rnaseq-custom.config run https://github.com/nf-core/rnaseq \
			-params-file $(CONFIG_DIR)/rnaseq-params.json \
			-profile test \
			--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml"

.PHONY: test_rnaseq_full
test_rnaseq_full: setup_environment
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(BASEDIR) && \
		spack env activate -d $(ENV_DIR) && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(CONFIG_DIR)/rnaseq-params.json \
			-profile test_full \
			--multiqc_config $(CONFIG_DIR)/custom_multiqc_config.yml"

.PHONY: clean
clean:
	@echo "Cleaning up..."
	@rm -rf spack ~/.spack
	@if [ -d "$(ENV_DIR)/.spack-env" ]; then \
		rm -rf "$(ENV_DIR)/.spack-env"; \
	fi
	@if [ -f "$(ENV_DIR)/spack.lock" ]; then \
		rm -f "$(ENV_DIR)/spack.lock"; \
	fi
