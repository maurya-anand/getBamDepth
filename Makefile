# Detect conda/mamba: prefer mamba, fall back to conda
CONDA := $(shell command -v mamba >/dev/null 2>&1 && echo mamba || (command -v conda >/dev/null 2>&1 && echo conda || echo ""))
ENV_NAME := getbamdepth

install:
	@if [ -z "$(CONDA)" ]; then \
		echo "ERROR: Neither mamba nor conda found in PATH."; \
		exit 1; \
	fi
	@chmod a+x ./getBamDepth
	@echo "Installing samtools and perl modules via conda..."
	@echo "Using $(CONDA) to create $(ENV_NAME) environment..."
	@$(CONDA) create -y -n $(ENV_NAME) -c conda-forge -c bioconda samtools perl-sys-cpu
	@echo ""
	@echo "Installation complete!"
	@echo "Activate conda environment: conda activate $(ENV_NAME)"

uninstall:
	@if [ -z "$(CONDA)" ]; then \
		echo "ERROR: Neither mamba nor conda found in PATH."; \
		exit 1; \
	fi
	@echo "Removing conda environment $(ENV_NAME)..."
	@$(CONDA) env remove -n $(ENV_NAME) -y
	@echo "Finished uninstalling conda environment..."

.PHONY: tests
tests:
	@if [ -z "$(CONDA)" ]; then \
		echo "ERROR: Neither mamba nor conda found in PATH."; \
		echo "Install mamba/conda before running tests."; \
		exit 1; \
	fi
	@echo "Running test suite..."
	@echo ""
	@echo "Test 1: Depth file input with default thresholds..."
	@$(CONDA) run -n $(ENV_NAME) ./getBamDepth --bed example/example-targets.bed --depth example/sample.depth --output /tmp/test1.txt
	@diff /tmp/test1.txt example/sample.coverage.out.txt && echo "[PASS] Test 1 passed" || (echo "[FAIL] Test 1 failed"; exit 1)
	@echo ""
	@echo "Test 2: CRAM input with default thresholds..."
	@$(CONDA) run -n $(ENV_NAME) ./getBamDepth --bed example/example-targets.bed --bam example/sample.cram --output /tmp/test2.txt
	@diff /tmp/test2.txt example/sample.coverage.out.txt && echo "[PASS] Test 2 passed" || (echo "[FAIL] Test 2 failed"; exit 1)
	@echo ""
	@echo "Test 3: BAM input with custom thresholds (5x,10x)..."
	@$(CONDA) run -n $(ENV_NAME) ./getBamDepth --bed example/example-targets.bed --bam example/sample.bam --thresholds 5,10 --output /tmp/test3.txt
	@if grep -q "5x" /tmp/test3.txt && grep -q "10x" /tmp/test3.txt; then echo "[PASS] Test 3 passed"; else echo "[FAIL] Test 3 failed"; exit 1; fi
	@echo ""
	@echo "All tests passed!"
