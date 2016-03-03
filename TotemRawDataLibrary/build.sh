#! /bin/bash

# Script used by Jenkins continuous integration server in order to build the
# project.

AGENT_USERNAME=`whoami`
AGENT_KERBEROS_KEYTAB="/etc/$AGENT_USERNAME.keytab"
AGENT_EXECUTORS=4
AGENT_BUILD_LOG="$WORKSPACE/build.log"

set -e # Exit immediately if a simple command exits with a non-zero status.
set -x # Log all commands to stdout.
set -o pipefail # Return value of a pipeline as the value of the last command to
                # exit with a non-zero status, or zero if all commands in the
                # pipeline exit successfully.

# Shows directory details.
function show_current_directory_details() {
	echo "Content of directory '`pwd`':"
	ls -al
}

# Initializes Kerberos keytab for agent user.
function initialize_kerberos_keytab() {
	echo "Initializing Kerberos keytab for user '$AGENT_USERNAME'..."
	kinit -kt "$AGENT_KERBEROS_KEYTAB" "$AGENT_USERNAME"
}

# Compiles cmake project using AGENT_EXECUTORS number of threads.
# $1 - CMake options
function compile_cmake_project() {
	echo "Compiling cmake project..."
	show_current_directory_details
	mkdir build
	cd build
	cmake .. $1
	show_current_directory_details
	make -j"$AGENT_EXECUTORS"
}

# Builds TotemRawDataLibrary project.
# $1 - URL to project SVN repository
function build_TotemRawDataLibrary() {
	echo "Building TotemRawDataLibrary project..."
	initialize_kerberos_keytab
	compile_cmake_project
}

# Executes given command and logs stdout both to the file and the screen.
# $1 - commnad to be executed
function execute_and_log() {
	echo "Executing command '$1'"
	eval "$1" 2>&1 | tee -a "$AGENT_BUILD_LOG"
}

execute_and_log "build_TotemRawDataLibrary"
