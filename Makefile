# Detect conda/mamba: prefer mamba, fall back to conda
CONDA := $(shell command -v mamba >/dev/null 2>&1 && echo mamba || (command -v conda >/dev/null 2>&1 && echo conda || echo ""))
ENV_NAME := bd-env
TEST_DIR := /tmp/bd-tests

.PHONY: help install uninstall test test-full

help:
	@echo "getBamDepth build targets:"
	@echo "  make install    - Install dependencies in conda environment"
	@echo "  make uninstall  - Remove conda environment"
	@echo "  make test       - Quick smoke tests (3 functional tests)"
	@echo "  make test-full  - Comprehensive tests (functional + validation)"

install:
	@if [ -z "$(CONDA)" ]; then \
		echo "ERROR: Neither mamba nor conda found in PATH."; \
		exit 1; \
	fi
	@chmod a+x ./getBamDepth
	@echo "Installing samtools and perl modules via conda..."
	@echo "Using $(CONDA) to create $(ENV_NAME) environment..."
	@if $(CONDA) info --envs | grep -q "$(ENV_NAME)"; then \
		echo "Environment $(ENV_NAME) already exists. Skipping creation."; \
	else \
		$(CONDA) create -y -n $(ENV_NAME) -c conda-forge -c bioconda samtools perl; \
	fi
	@echo "Installing getBamDepth to environment bin directory..."
	@CONDA_BIN=$$($(CONDA) run -n $(ENV_NAME) sh -c 'echo $$CONDA_PREFIX'/bin); \
	mkdir -p "$$CONDA_BIN" 2>/dev/null || true; \
	cp ./getBamDepth "$$CONDA_BIN/getBamDepth"; \
	chmod +x "$$CONDA_BIN/getBamDepth"
	@echo ""
	@echo "Installation complete!"
	@echo "Activate conda environment: conda activate $(ENV_NAME)"
	@echo "Then run: getBamDepth --help"

uninstall:
	@if [ -z "$(CONDA)" ]; then \
		echo "ERROR: Neither mamba nor conda found in PATH."; \
		exit 1; \
	fi
	@echo "Removing conda environment $(ENV_NAME)..."
	@$(CONDA) env remove -n $(ENV_NAME) -y
	@echo "Finished uninstalling conda environment..."

test:
	@if [ -z "$(CONDA)" ]; then \
		echo "ERROR: Neither mamba nor conda found in PATH."; \
		echo "Install mamba/conda before running tests."; \
		exit 1; \
	fi
	@mkdir -p $(TEST_DIR); \
	set -e; \
	trap "rm -rf $(TEST_DIR)" EXIT; \
	echo "Running quick smoke tests..."; \
	echo ""; \
	echo "Test 1: Depth file input with default thresholds..."; \
	$(CONDA) run -n $(ENV_NAME) ./getBamDepth --bed example/example-targets.bed --depth example/sample.depth --output $(TEST_DIR)/test1.txt && \
	diff $(TEST_DIR)/test1.txt example/sample.coverage.out.txt && echo "[PASS] Test 1 passed" || { echo "[FAIL] Test 1 failed"; exit 1; }; \
	echo ""; \
	echo "Test 2: CRAM input with default thresholds..."; \
	$(CONDA) run -n $(ENV_NAME) ./getBamDepth --bed example/example-targets.bed --bam example/sample.cram --output $(TEST_DIR)/test2.txt && \
	diff $(TEST_DIR)/test2.txt example/sample.coverage.out.txt && echo "[PASS] Test 2 passed" || { echo "[FAIL] Test 2 failed"; exit 1; }; \
	echo ""; \
	echo "Test 3: BAM input with custom thresholds (5x,10x)..."; \
	$(CONDA) run -n $(ENV_NAME) ./getBamDepth --bed example/example-targets.bed --bam example/sample.bam --thresholds 5,10 --output $(TEST_DIR)/test3.txt && \
	grep -q "5x" $(TEST_DIR)/test3.txt && grep -q "10x" $(TEST_DIR)/test3.txt && echo "[PASS] Test 3 passed" || { echo "[FAIL] Test 3 failed"; exit 1; }; \
	echo ""; \
	echo "All smoke tests passed!"

test-full: test
	@echo ""
	@bash tests/validation.sh "$(CONDA)" "$(ENV_NAME)"
