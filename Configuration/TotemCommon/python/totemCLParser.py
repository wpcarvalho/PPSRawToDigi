"""
 Totem command line parser
 this parser defines several common options such as nevents, beta,verbosity, help, energy etc.
 for documentation see also: https://twiki.cern.ch/twiki/bin/view/TOTEM/TotemCommandLineParser
"""
import sys as _sys

# for details of python argparse modul see: http://code.google.com/p/argparse/
import Configuration.TotemCommon.argparse as _argparse

class TotemCLParser(_argparse.ArgumentParser):
  #def __init__(self,*args, **kw):
  def __init__(self):
    # this if/else switch allows to use this parser in both of the following situations
    if "cmsRun" in _sys.argv:
      progName = _sys.argv[1] # if a script is run like: cmsRun script.py +options (_sys.argv[1]=="cmsRun")
    else:
      progName = _sys.argv[0] # if a script is run like: python script.py +options (_sys.argv[0]=="script.py")

    _argparse.ArgumentParser.__init__(self,
        # prog=_sys.argv[1],  # default name of "program" (_sys.argv[0]="cmsRun")
        prog=progName,      # default name of program
        prefix_chars='-+',  # -options with dash ('-') prefix are recognized by cmsRun and options with '+' prefix are recognized by a configuration file 
        add_help=False      # do not add default help (it contains -h --help options which are handled by cmsRun and not by this parser)
        )

    # this positional argument has to be defined to suppres one error... (if cmsRun is used)
    if "cmsRun" in _sys.argv:
      self.add_argument(
        "cfgfilename", 
        default=self.prog, 
        action='store', 
        help=_argparse.SUPPRESS #help = ""
      )

    self.add_argument(
        "+nevents",           # option name(s) used from command line
        dest    = 'nevents',  # one common name
        default = None,       # dafault value
        #default = -1,        # dafault value
        #    default = _argparse.SUPPRESS,       # dafault value
        type    = int,      
        action  = "store",    
        help    = "number of events to be processed"
        )

    self.add_argument(
        "+beta",
        dest    = 'beta',
        default = None,
        action  = "store",
        help    = "beta* - optics"
        )

    self.add_argument(
        "+energy",
        dest    = 'energy',
        default = None,
        type    = float,
        action  = "store",
        help    = "energy"
        )

    self.add_argument(
        "+verbosity",
        dest    = 'verbosity',
        default = None,
        type    = int,
        action  = "store",
        help    = "verbosity"
        )

    self.add_argument(
        '+h', '++help', 
        action  = 'help',
        default = _argparse.SUPPRESS,
        help    = 'show this help message and exit'
        )

    #self.add_argument(
#   '+v', '++version', 
#   action  = 'version', 
#   default = _argparse.SUPPRESS,
#   version = _argparse.version,
#   help    = "show program's version number and exit"
#)
  def parse_args(self,*args, **kw):
    # redefine this function so that parsed arguments are also printed to the screen
    args = _argparse.ArgumentParser.parse_args(self,*args, **kw)
    print "Parsed CL arguments: \n", args, "\n"
    return  args

parser = TotemCLParser()




