#!/bin/bash
#
# Script for running code coverage of ldpred tests.
#
# For this script to work, the coverage tool must have previously been installed
# via `pip install coverage` .
#
# Coverage files will be stored in subdirectories under ${TMPDIR}/ldpred or, if
# TMPDIR is not defined, /tmp/ldpred.
#
# By default, all tests are run.  However, any arguments passed to coverage.sh
# will be forwarded to the `coverage` tool.  See below for an example on how to
# run a specific test.
#
# Examples
# --------
# Run coverage on all tests, using default tmp dir.  This will create /tmp/ldpred/.
# $ ./coverage.sh
#
# Run coverage on all tests , using specified tmp dir.  This will create ${HOME}/tmp/ldpred.
# $ TMPDIR=${HOME}/tmp ./coverage.sh
#
# Run coverage for a specific test.
# $ ./coverage.sh -m unittest test.TestLDpred.test_ldpred_inf

# Check first if the coverage tool is available.
readonly COVERAGE_BIN="$(which coverage)"
if [[ ! -x "${COVERAGE_BIN}" ]]; then
  echo
  echo "**** The 'coverage' tool is not available on the path."
  echo "**** To install coverage, run:"
  echo "****     pip install coverage"
  echo "****"
  echo "**** For details, see https://coverage.readthedocs.io"
  echo
  exit 1
fi

# Get fully qualified path of the ldpred source directory.
readonly LDPRED_SRC_ROOT="$(readlink -e "$(dirname "$0")")"
echo "Running code coverage for ldpred sources under: ${LDPRED_SRC_ROOT}"

# Create temporary directory for coverage output.
readonly COVERAGE_TMPDIR="${TMPDIR:-/tmp}/ldpred/coverage_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${COVERAGE_TMPDIR}"
cd "${COVERAGE_TMPDIR}"
echo "The coverage temp directory is: ${COVERAGE_TMPDIR}"

# Run code coverage and generate HTML report.
PYTHONPATH="${LDPRED_SRC_ROOT}:${PYTHONPATH}" python -B "${COVERAGE_BIN}" run --branch --source="${LDPRED_SRC_ROOT}" "${@:-${LDPRED_SRC_ROOT}/test.py}"
coverage html
echo "To view the code coverage report, point your web browser to:"
echo "file:${COVERAGE_TMPDIR}/htmlcov/index.html"
echo
