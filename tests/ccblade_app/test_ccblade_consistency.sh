#!/bin/bash

# Script to prepare test data for CCBlade consistency tests

TEST_DIR="tests/ccblade_app/test_data"
OUTPUT_DIR="$TEST_DIR/output"
REFERENCE_DIR="$TEST_DIR/reference"
TEST_YAML="$TEST_DIR/test_config.yml"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REFERENCE_DIR"

# Check if test config exists
if [ ! -f "$TEST_YAML" ]; then
    echo "Test config $TEST_YAML not found. Please provide it."
    exit 1
fi

# Run the CCBlade app to generate outputs
b3p ccblade -y "$TEST_YAML" run

# Copy outputs to reference directory
if [ -d "$OUTPUT_DIR" ] && [ "$(ls -A $OUTPUT_DIR)" ]; then
    cp "$OUTPUT_DIR"/* "$REFERENCE_DIR"/
    echo "Reference data prepared in $REFERENCE_DIR"
else
    echo "No outputs found in $OUTPUT_DIR"
fi
