TOTEM Raw Data Package
====================



Access
---------------------

The source code can be accessed via TOTEM Subversion repository (svn+ssh://svn.cern.ch/reps/totem/trunk/TotemRawData).
You can view code using [WebSVN](https://svnweb.cern.ch/cern/wsvn/totem/trunk/TotemRawData)



Technologies
---------------------

+ Castor (rfio)
+ CMake
+ optionally GMock



Supported versions
---------------------

+ **gcc v4.8.1+** as long you have libraries built with older cpp standard library you may use older gcc versions; (if you need it on lxplus you can setup cmssw v7.0.4+)
+ **Castor v2.1.13-6** only that version was tested; please note that rfio obsolete and is required by project for the time being; if no castor is found monitor will adapt and buid without Castor support
+ **CMake v2.8+** some macros introduced here are needed; on lxplus (with SLC6) use command **cmake28**


Output
---------------------

Produces library files which you can link to your application.
As well some testing binaries are produced (which you can run manually).
As an output we also consider GTest test cases (which you can launch manually or via ctest).



Usage
---------------------

To build project use standard cmake procedure

```sh
mkdir <workingdir> && cd <workingdir>
cmake ..
make
```

If you are interested in using only monitor binary use make target monitor.
Please take a note that CMake may not find proper library versions. You may consider specifying them manually by passing definitions to cmake

```sh
cmake .. \
-DCASTOR_DIR=<Castor build dir> \
[-GMOCK_ROOT=<GMock library build root (includes GTest)> (without that tests won't be available)] \
```

which could look like

```sh
cmake .. \
-DCASTOR_DIR=/afs/cern.ch/sw/lcg/external/castor/2.1.13-6/x86_64-slc6-gcc47-opt/usr/
-GMOCK_ROOT=(not on afs)
```

alternatively you can choose libraries versions by setting environment variables

```sh
export CASTOR_DIR="/afs/cern.ch/sw/lcg/external/castor/2.1.13-6/x86_64-slc6-gcc47-opt/usr/"
```
