
fold:
	cfold fold 

foldt1:
	cfold fold b3p/**py tests/*py


unfold:
	cfold unfold v1.txt
	pip install -e . 






.PHONY: test test-verbose test-coverage clean


# Test directories and files
TEST_FILES = tests/test_2d.py  tests/test_bondline_model.py  tests/test_build.py  tests/test_drape.py  tests/test_geometry.py tests/test_ccx.py

# Default target
all: test

# Run tests
test:
	poetry run pytest $(TEST_FILES)

# Run tests with verbose output
test-verbose:
	poetry run pytest $(TEST_FILES) -v

# Run tests with coverage report
test-coverage:
	poetry run pytest --cov=b3p --cov-report=term --cov-report=html


# Clean up generated files
clean:
	rm -rf __pycache__
	rm -rf .pytest_cache
	rm -rf htmlcov
	rm -rf .coverage
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete
