#! /bin/bash

# Script used by Jenkins continuous integration server in order to build the
# project.

AGENT_EXECUTORS=4
AGENT_BUILD_LOG="$WORKSPACE/build.log"

set -e # Exit immediately if a simple command exits with a non-zero status.
set -x # Log all commands to stdout.
set -o pipefail # Return value of a pipeline as the value of the last command to
                # exit with a non-zero status, or zero if all commands in the
                # pipeline exit successfully.
shopt -s expand_aliases # Expand command alias to the command itself.
                        # Required for non-interactive shell.

if [ -z "$SCRAM_ARCH" ]; then
    export SCRAM_ARCH=slc6_amd64_gcc493
fi
source /afs/cern.ch/cms/cmsset_default.sh

# Shows directory details.
function show_current_directory_details() {
	echo "Content of directory '`pwd`':"
	ls -al
}

# Initializes scram project.
# $1 - project name
# $2 - project version
function initialize_scram_project() {
	echo "Initializing scram project..."
	show_current_directory_details
	scram project -n "$JOB_NAME" "$1" "$2"
	show_current_directory_details
	cd "$JOB_NAME"
	show_current_directory_details
	cmsenv
}

# Compiles scram project using AGENT_EXECUTORS number of threads.
function compile_scram_project() {
	echo "Compiling scram project..."
	cd src
	show_current_directory_details
	echo "Starting parallel compilation using up to $AGENT_EXECUTORS threads..."
	scram build -j"$AGENT_EXECUTORS" || \
	scram build -j"$AGENT_EXECUTORS" || \
	echo "Parallel compilation using up to $AGENT_EXECUTORS threads failed..."
	echo "Starting sequential compilation..."
	scram build
}

# Builds CMSSW project, version 8.0.0.
# $1 - if argument equals "buildMergingSoftware" an additional build step will
#      be executed
function build_CMSSW_8_0_0() {
	echo "Building CMSSW 8.0.0 project..."
	initialize_scram_project "CMSSW" "CMSSW_8_0_0_pre5"
	compile_scram_project
}

# Executes given command and logs stdout both to the file and the screen.
# $1 - commnad to be executed
function execute_and_log() {
	echo "Executing command '$1'"
	eval "$1" 2>&1 | tee -a "$AGENT_BUILD_LOG"
}

execute_and_log "build_CMSSW_8_0_0 $1"
