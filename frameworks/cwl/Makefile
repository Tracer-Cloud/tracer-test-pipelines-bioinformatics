POETRY := /root/.local/bin/poetry

.PHONY: install_poetry setup_env create_inputs run_cwl

install_poetry:
	@if ! command -v $(POETRY) &> /dev/null; then \
		echo "Installing Poetry..."; \
		curl -sSL https://install.python-poetry.org | python3 -; \
	else \
		echo "Poetry already installed."; \
	fi

setup_env:
	$(POETRY) install

create_inputs:
	bash generate.sh

run_cwl: install_poetry setup_env create_inputs
	$(POETRY) run cwltool rnaseq-workflow.cwl rnaseq-job-001.yml
