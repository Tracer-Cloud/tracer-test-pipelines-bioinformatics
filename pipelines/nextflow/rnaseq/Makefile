SHELL := bash
SPACK_DIR := $(abspath $(CURDIR)/spack)
ENV_DIR := $(abspath $(CURDIR)/..)
CONFIG_DIR := $(abspath $(ENV_DIR)/config)
PARAMS_DIR := $(abspath $(CONFIG_DIR))
MULTIQC_CONFIG := $(abspath $(CONFIG_DIR)/custom_multiqc_config.yml)
TEMP_PARAMS := $(abspath $(CURDIR)/temp_params.json)

# Get parent of makefile which 
BASEDIR := $(abspath $(CURDIR)/)

# Function to create dynamic params file
define create_dynamic_params
	@echo "Creating dynamic params file..."
	@echo '{' > $(TEMP_PARAMS)
	@echo '    "outdir": "$(1)",' >> $(TEMP_PARAMS)
	@echo '    "multiqc_config": "$(MULTIQC_CONFIG)"' >> $(TEMP_PARAMS)
	@echo '}' >> $(TEMP_PARAMS)
endef

# Function to clean up temp params file
define cleanup_temp_params
	@rm -f $(TEMP_PARAMS)
endef

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
	@if ! grep -q "source.*spack/setup-env.sh" ~/.bashrc; then \
		echo "source ~/.spack/setup-env.sh" >> ~/.bashrc; \
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
		if [ ! -d "$(ENV_DIR)/.spack-env" ]; then \
			cd "$(ENV_DIR)" && \
			spack env create -d . --with-view "$(ENV_DIR)/.spack-env/view" && \
			spack env activate -d . && \
			spack install; \
		fi'
	@echo "Spack environment created successfully."

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
	$(call create_dynamic_params,output/rnaseq-test)
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(TEMP_PARAMS) \
			-profile test"
	$(call cleanup_temp_params)

.PHONY: test_rnaseq_full
test_rnaseq_full: setup_environment
	$(call create_dynamic_params,output/rnaseq-test)
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(TEMP_PARAMS) \
			-profile test_full"
	$(call cleanup_temp_params)

.PHONY: test_rnaseq_right_size_default
test_rnaseq_right_size_default: setup_environment
	$(call create_dynamic_params,s3://tracer-nxf-outputs/rnaseq/right-size-tests/default)
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(TEMP_PARAMS) \
			-profile test_full"
	$(call cleanup_temp_params)

.PHONY: test_rnaseq_right_size_half_cpu
test_rnaseq_right_size_half_cpu: setup_environment
	$(call create_dynamic_params,s3://tracer-nxf-outputs/rnaseq/right-size-tests/half-cpu)
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(TEMP_PARAMS) \
			-profile test_full"
	$(call cleanup_temp_params)

.PHONY: test_rnaseq_right_size_double_cpu
test_rnaseq_right_size_double_cpu: setup_environment
	$(call create_dynamic_params,s3://tracer-nxf-outputs/rnaseq/right-size-tests/double-cpu)
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
		cd $(ENV_DIR) && \
		spack env activate -d . && \
		nextflow -c $(CONFIG_DIR)/local.config run https://github.com/nf-core/rnaseq \
			-params-file $(TEMP_PARAMS) \
			-profile test_full"
	$(call cleanup_temp_params)

.PHONY: clean
clean:
	@echo "Cleaning up..."
	@rm -rf spack ~/.spack $(TEMP_PARAMS)
