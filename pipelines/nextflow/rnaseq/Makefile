SHELL := bash
SPACK_DIR := $(abspath $(CURDIR)/spack)
ENV_DIR := $(abspath $(CURDIR)/..)
CONFIG_DIR := $(abspath $(ENV_DIR)/config)

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
		echo "Cloning Spack repository (this may take a few minutes)..."; \
		git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git > /dev/null 2>&1; \
		echo "Spack repository cloned successfully."; \
	fi

.PHONY: setup_spack_env
setup_spack_env: install_spack
	@echo "Setting up Spack environment..."
	@if [ ! -f /etc/profile.d/spack.sh ]; then \
		if [ "$$(id -u)" = "0" ]; then \
			echo "Creating system-wide Spack configuration..."; \
			echo "export SPACK_ROOT=$(SPACK_DIR)" > /etc/profile.d/spack.sh; \
			echo "source \$$SPACK_ROOT/share/spack/setup-env.sh" >> /etc/profile.d/spack.sh; \
			echo "export PATH=\$$SPACK_ROOT/bin:\$$PATH" >> /etc/profile.d/spack.sh; \
			chmod +x /etc/profile.d/spack.sh; \
		else \
			echo "Creating system-wide Spack configuration (requires sudo)..."; \
			echo "export SPACK_ROOT=$(SPACK_DIR)" | sudo tee /etc/profile.d/spack.sh > /dev/null; \
			echo "source \$$SPACK_ROOT/share/spack/setup-env.sh" | sudo tee -a /etc/profile.d/spack.sh > /dev/null; \
			echo "export PATH=\$$SPACK_ROOT/bin:\$$PATH" | sudo tee -a /etc/profile.d/spack.sh > /dev/null; \
			sudo chmod +x /etc/profile.d/spack.sh > /dev/null; \
		fi \
	fi
	@echo "Configuring user environment..."
	@if ! grep -q "source.*spack/setup-env.sh" ~/.bashrc; then \
		echo "export SPACK_ROOT=$(SPACK_DIR)" >> ~/.bashrc; \
		echo "source \$$SPACK_ROOT/share/spack/setup-env.sh" >> ~/.bashrc; \
		echo "export PATH=\$$SPACK_ROOT/bin:\$$PATH" >> ~/.bashrc; \
		source ~/.bashrc > /dev/null 2>&1; \
	fi
	@echo "Verifying Spack installation..."
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && spack --version" > /dev/null 2>&1 || { \
		echo "Error: Failed to set up Spack environment"; \
		exit 1; \
	}
	@echo "Spack environment setup completed successfully."

.PHONY: verify_spack
verify_spack: setup_spack_env
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		spack --version > /dev/null 2>&1 && \
		spack help --all > /dev/null 2>&1" || { \
		echo "Error: Spack verification failed"; \
		exit 1; \
	}

.PHONY: create_env
create_env: verify_spack
	@echo "Creating Spack environment..."
	@bash -c 'source "$(SPACK_DIR)/share/spack/setup-env.sh" && \
		if [ ! -d "$(ENV_DIR)/.spack-env" ]; then \
			cd "$(ENV_DIR)" && \
			echo "Using existing spack.yaml from root directory..."; \
			spack env create -d . --with-view "$(ENV_DIR)/.spack-env/view" && \
			echo "Activating environment..."; \
			spack env activate -d . && \
			echo "Installing packages from spack.yaml..."; \
			spack install; \
		fi'
	@echo "Spack environment created successfully."

.PHONY: setup_environment
setup_environment: create_env
	@echo "Installing dependencies (this may take several minutes)..."
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . > /dev/null 2>&1 && \
		spack install > /dev/null 2>&1"
	@echo "Environment setup completed successfully."

.PHONY: test_rnaseq
test_rnaseq: setup_environment
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(CONFIG_DIR)/rnaseq-params.json \
			-profile test"

.PHONY: test_rnaseq_full
test_rnaseq_full: setup_environment
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(CONFIG_DIR)/rnaseq-params.json \
			-profile test_full"

.PHONY: clean
clean:
	@echo "Cleaning up..."
	@rm -rf spack
	@if [ -f /etc/profile.d/spack.sh ]; then \
		sudo rm /etc/profile.d/spack.sh; \
	fi
