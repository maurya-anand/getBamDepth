#!/bin/bash

# Comprehensive validation test suite for getBamDepth

CONDA="${1:-conda}"
ENV_NAME="${2:-getbamdepth}"
SCRIPT="./getBamDepth"
DATA_DIR="tests/data"
OUTPUT_DIR="/tmp/getBamDepth-validation-tests"
PASS=0
FAIL=0

# Prepare command wrapper to run in conda environment
CMD_PREFIX="$CONDA run -n $ENV_NAME"

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

cleanup() {
    [ -d "$OUTPUT_DIR" ] && rm -rf "$OUTPUT_DIR" 2>/dev/null || true
}
trap cleanup EXIT

mkdir -p "$OUTPUT_DIR" 2>/dev/null || true

test_pass() {
    echo -e "${GREEN}✓${NC} $1"
    ((PASS++))
}

test_fail() {
    echo -e "${RED}✗${NC} $1"
    ((FAIL++))
}

run_test() {
    local name="$1" cmd="$2" should_fail="$3" pattern="$4"
    local output exit_code

    output=$(bash -c "$cmd" 2>&1)
    exit_code=$?

    if [ "$should_fail" = "true" ]; then
        if [ $exit_code -ne 0 ]; then
            if [ -n "$pattern" ]; then
                if echo "$output" | grep -q "$pattern"; then
                    test_pass "$name"
                else
                    test_fail "$name (wrong error)"
                fi
            else
                test_pass "$name"
            fi
        else
            test_fail "$name (should fail)"
        fi
    else
        if [ $exit_code -eq 0 ]; then
            test_pass "$name"
        else
            test_fail "$name"
        fi
    fi
}

echo "=========================================="
echo "getBamDepth Validation Test Suite"
echo "=========================================="
echo ""

echo ">>> THRESHOLD VALIDATION <<<"
run_test "Accept valid integers" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --thresholds 5,10 --output $OUTPUT_DIR/out1.txt" false ""
run_test "Accept valid floats" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --thresholds 5.5,10.5 --output $OUTPUT_DIR/out2.txt" false ""
run_test "Reject non-numeric" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --thresholds abc,def --output $OUTPUT_DIR/out3.txt" true "Invalid threshold"
run_test "Reject negative" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --thresholds -5,10 --output $OUTPUT_DIR/out4.txt" true "Invalid threshold"
run_test "Reject zero" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --thresholds 0,10 --output $OUTPUT_DIR/out5.txt" true "Invalid threshold"
run_test "Accept unsorted" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --thresholds 20,5,10 --output $OUTPUT_DIR/out6.txt" false ""
echo ""

echo ">>> BED FILE VALIDATION <<<"
run_test "Accept valid BED" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --output $OUTPUT_DIR/out7.txt" false ""
run_test "Reject <3 columns" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/bad_cols.bed --depth $DATA_DIR/good.depth --output $OUTPUT_DIR/out8.txt" true "expected at least 3 columns"
run_test "Reject non-integer coords" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/bad_coords.bed --depth $DATA_DIR/good.depth --output $OUTPUT_DIR/out9.txt" true "non-negative integer"
run_test "Reject inverted range" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/inverted.bed --depth $DATA_DIR/good.depth --output $OUTPUT_DIR/out10.txt" true "end must be >= start"
echo ""

echo ">>> DEPTH FILE VALIDATION <<<"
run_test "Accept valid depth" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/good.depth --output $OUTPUT_DIR/out11.txt" false ""
run_test "Reject <3 columns" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/bad_depth_cols.depth --output $OUTPUT_DIR/out12.txt" true "expected 3 columns"
run_test "Reject non-integer pos" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/bad_pos.depth --output $OUTPUT_DIR/out13.txt" true "positive integer"
run_test "Reject non-integer depth" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/good.bed --depth $DATA_DIR/bad_depth_val.depth --output $OUTPUT_DIR/out14.txt" true "non-negative integer"
echo ""

echo ">>> EDGE CASES <<<"
run_test "Handle empty region" "$CMD_PREFIX $SCRIPT --bed $DATA_DIR/empty.bed --depth $DATA_DIR/good.depth --output $OUTPUT_DIR/out15.txt" false ""
echo ""

echo "=========================================="
echo "VALIDATION SUMMARY"
echo "=========================================="
echo -e "Passed: ${GREEN}$PASS${NC}"
echo -e "Failed: ${RED}$FAIL${NC}"
echo ""

if [ $FAIL -eq 0 ]; then
    echo -e "${GREEN}All 16 validation tests passed! ✓${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed.${NC}"
    exit 1
fi
