RPStationMultiTrackFinderFitter
===============================


Configurations
-------------
Simulation is to be runned from the test directory. There is a stack of
configurations that are used.

Geometry, track distribution are defined in configurations (usually named
`simulate_cfg.py`) placed in test/**scenario_name** directory. Each
config file can be seen as an ordinary Python module.

There are also modules at test directory level which do not modify distribution
or reconstruction parameters but adjust the way simulations are run and what
output is generated.

They operate on a track distribution which is passed to them as a name of
a Python module.

### Statistical test of performance of a large number of events.
```
cmsRun set_particle_count_cfg.py <module> <track count> <event count> [<RP> ...]
```
For example, in order to run a 4000GeV Beta\*=90m track distribution with two
tracks present and 1000 events generated using only RPs with IDs: 100, 104
and 120 one should run
```
cmsRun set_particle_count_cfg.py 4000GeV_90m.simulate_cfg 2 1000 100 104 120
```

Verbose output is suppressed.
Output is generated in stats_<track count> and hists_<track count>.

stats\_<track count> has following format:
```
<run>:<event> <missing>/<fake>
<run>:<event> <missing>/<fake>
...
```
6 last lines contain statistical analysis of the former lines, they should
be self explanatory.

### A thorough analysis of a single event
```
cmsRun selected_event_analysis_cfg.py <module> <track count> <event number> [<RP> ...]
```
No statistics are generated in such run. All verbose output is enabled what
makes it possible to analyse what went wrong (or right) in selected event.

To understand the output reading the docs and the code is advised.

Should the need arise it is also possible to pipe the output to 
`utils/verbose_to_asy.py [RPid ...]` in order to obtain visualisation of event.
The output of the tool in in asymptote code so it should be rendered with
`asy` in order to obtain readable graphs.

### Generate track distribution
There is also a `generate_cfg.py` configuration used to generate the track 
distributions as it may be more time consuming than the reconstruction 
itself (details on how to use it can be found in each distribution directory).

Documentation
-------------
In order to generate code documentation run (in module root directory):
```
doxygen Doxyfile
```
output will be written to `doc` directory.
