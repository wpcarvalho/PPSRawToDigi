TOTEM Raw Data Package
====================



Access
---------------------

The source code can be accessed via TOTEM Subversion repository (svn+ssh://svn.cern.ch/reps/totem/trunk/TotemRawData).
You can view code using [WebSVN](https://svnweb.cern.ch/cern/wsvn/totem/trunk/TotemRawData)



Technologies
---------------------
+ c++11
+ ROOT
+ CMake
+ Castor (rfio)
+ CMake



Supported versions
---------------------
+ **gcc v4.8.1+** required libraries like ROOT builds are performed with such; as long you have libraries built with older cpp standard library you may use older gcc versions; (if you need it on lxplus you can setup cmssw v7.0.4+)
+ **ROOT v5.34.x** versions up to 5.34.19 were tested; there is no guarantee that other versions would work
+ **Castor v2.1.13-6** only that version was tested; please note that rfio obsolete and is required by project for the time beeing; if no castor is found monitor will adapt and buid without Castor support
+ **CMake v2.8+** some macros introduced here are needed; on lxplus (with SLC6) use command **cmake28**


Output
---------------------

Produces library files which you can link to your application.
As well some testing binaries are produced (which you can run manuallly).
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
Please take a note that CMake may not find proper library versions. You may consider specyfying them manually by passing definitions to cmake

```sh
cmake .. \
-DQT_QMAKE_EXECUTABLE=<qmake binary path> \
-DLibXML2_ROOT_DIR=<libxml2 build dir> \
-DROOT_CONFIG_SEARCHPATH=<ROOT binary dir> \
-DCASTOR_DIR=<Castor build dir> \
[-GMOCK_ROOT=<GMock library build root (includes GTest)> (without that tests won't be available)] \
```

which could look like

```sh
cmake .. \
-DQT_QMAKE_EXECUTABLE=/afs/cern.ch/sw/lcg/external/qt/4.8.4/x86_64-slc6-gcc47-opt/bin/qmake \
-DLibXML2_ROOT_DIR=/afs/cern.ch/exp/totem/soft/libxml/libxml2-2.9.0-x86_64-gcc4.4 \
-DROOT_CONFIG_SEARCHPATH=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.17/x86_64-slc6-gcc47-dbg/root/bin \
-DCASTOR_DIR=/afs/cern.ch/sw/lcg/external/castor/2.1.13-6/x86_64-slc6-gcc47-opt/usr/
-GMOCK_ROOT=(not on afs)
```

altirnatevely you can choose libraries versions by setting environment variables

```sh
export ROOTSYS="/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.17/x86_64-slc6-gcc47-dbg/root"
export QTDIR="/afs/cern.ch/sw/lcg/external/qt/4.8.4/x86_64-slc6-gcc47-opt"
export LIBXML2_DIR="/afs/cern.ch/exp/totem/soft/libxml/libxml2-2.9.0-x86_64-gcc4.4"
export CASTOR_DIR="/afs/cern.ch/sw/lcg/external/castor/2.1.13-6/x86_64-slc6-gcc47-opt/usr/"
```

notice that those variables point to slightly different locations - it is so because some of them have meaning in contexts different from monitor build (like ROOTSYS or QTDIR)