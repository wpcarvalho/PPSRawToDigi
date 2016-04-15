#!/usr/bin/env bash

set -e # Exit immediately if a simple command exits with a non-zero status.
set -x # Log all commands to stdout.
set -o pipefail # Return value of a pipeline as the value of the last command to
                # exit with a non-zero status, or zero if all commands in the
                # pipeline exit successfully.
shopt -s expand_aliases # Expand command alias to the command itself.
                        # Required for non-interactive shell.

source /opt/cms/cmsset_default.sh
scram b
eval `scram runtime -sh`
cmsRun raw_data_chain_test.py