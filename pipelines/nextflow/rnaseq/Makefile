SHELL := bash
SPACK_DIR := $(abspath $(CURDIR)/spack)
ENV_DIR := $(abspath $(CURDIR)/../../..)
CONFIG_DIR := $(abspath $(ENV_DIR)/pipelines/nextflow/config)
DEPS_SCRIPT := $(abspath $(ENV_DIR)/dependencies_installation.sh)

# Get parent of makefile which 
BASEDIR := $(abspath $(CURDIR)/)

.PHONY: check_prerequisites
check_prerequisites:
	@if [ ! -f "$(DEPS_SCRIPT)" ]; then \
		echo "Error: dependencies_installation.sh not found at $(DEPS_SCRIPT)"; \
		exit 1; \
	fi
	@bash "$(DEPS_SCRIPT)"

.PHONY: install_spack
install_spack: check_prerequisites
	@if [ ! -d "$(SPACK_DIR)" ]; then \
		echo "Installing Spack..."; \
		git clone --depth=1 https://github.com/spack/spack.git > /dev/null 2>&1; \
		echo "Spack repository cloned successfully."; \
	fi

.PHONY: setup_spack_env
setup_spack_env: install_spack
	@echo "Setting up Spack..."
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
	@echo "Spack setup completed successfully."

.PHONY: verify_spack
verify_spack: setup_spack_env
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && spack --version" > /dev/null 2>&1 || { \
		echo "Error: Spack verification failed"; \
		exit 1; \
	}

.PHONY: setup_compiler
setup_compiler: verify_spack
	@echo "Setting up compiler..."
	@bash -c 'source "$(SPACK_DIR)/share/spack/setup-env.sh" && \
		if command -v gcc &> /dev/null; then \
			echo "Using GCC compiler"; \
			spack compiler find --scope site; \
		elif command -v clang &> /dev/null; then \
			echo "Using Clang compiler"; \
			spack compiler find --scope site; \
		else \
			echo "Installing GCC compiler..."; \
			if command -v yum &> /dev/null; then \
				sudo yum install -y gcc gcc-c++ make; \
			elif command -v apt-get &> /dev/null; then \
				sudo apt-get install -y gcc g++ make; \
			fi; \
			spack compiler find --scope site; \
		fi'

.PHONY: setup_environment
setup_environment: verify_spack
	@echo "Setting up Spack environment..."
	@bash -c 'source "$(SPACK_DIR)/share/spack/setup-env.sh" && \
		cd "$(ENV_DIR)" && \
		spack env activate -d . && \
		spack concretize && \
		spack install'
	@echo "Environment setup completed successfully."

.PHONY: test_rnaseq
test_rnaseq:
	nextflow run nf-core/rnaseq -c $(CURDIR)/../../../custom.config -profile docker,test --outdir results -resume

.PHONY: test_rnaseq_full
test_rnaseq_full:
	nextflow run nf-core/rnaseq -c $(CURDIR)/../../../custom.config -profile docker,test --outdir results --full -resume

.PHONY: clean_spack
clean_spack:
	@echo "Cleaning Spack environment..."
	@bash -c 'source "$(SPACK_DIR)/share/spack/setup-env.sh" && \
		cd "$(ENV_DIR)" && \
		if [ -d ".spack-env" ]; then \
			rm -rf .spack-env; \
		fi && \
		if [ -f "spack.lock" ]; then \
			rm spack.lock; \
		fi'
	@echo "Spack environment cleaned successfully."

.PHONY: clean
clean: clean_spack
	@echo "Cleaning up..."
	@rm -rf spack ~/.spack

.PHONY: run_rnaseq_spack
run_rnaseq_spack:
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
	cd $(ENV_DIR) && \
	spack env activate -d . && \
	nextflow run nf-core/rnaseq -c $(CONFIG_DIR)/local.config -profile test --outdir results -resume"

.PHONY: run_rnaseq_spack_full
run_rnaseq_spack_full:
	@bash -c "source $(SPACK_DIR)/share/spack/setup-env.sh && \
	cd $(ENV_DIR) && \
	spack env activate -d . && \
	nextflow run nf-core/rnaseq -c $(CONFIG_DIR)/local.config -profile test --outdir results --full -resume"
