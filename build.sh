#!/usr/bin/env bash

set -e # Exit immediately if a simple command exits with a non-zero status.
set -x # Log all commands to stdout.
set -o pipefail # Return value of a pipeline as the value of the last command to
                # exit with a non-zero status, or zero if all commands in the
                # pipeline exit successfully.
shopt -s expand_aliases # Expand command alias to the command itself.
                        # Required for non-interactive shell.

cp -r src/* CMSSW_8_1_0_pre1/src
chown -R cmsbuild:cmsbuild CMSSW_8_1_0_pre1
source /opt/cms/cmsset_default.sh
cd CMSSW_8_1_0_pre1/src
scram b
eval `scram runtime -sh`
cmsRun raw_data_chain_test.py