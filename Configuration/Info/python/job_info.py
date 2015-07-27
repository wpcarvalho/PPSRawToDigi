# author: Leszek Grzanka (Leszek.Grzanka@cern.ch)

import sys
import os
import string
import getpass
import subprocess


# function to get actual username
def getUsername():
    result = getpass.getuser()
    return result


# function to get svn revision at given directory
def getSVNrevision( location ):
    program = ['svnversion']
    p = subprocess.Popen(program, shell=True, stdout=subprocess.PIPE, cwd=location, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    result = stdout.strip()
    if stderr.strip() != "":
        result += " (ERROR: " + stderr.strip() + ")"
    return result


# prints usage	
def usage():
    print "Usage:  %s" % os.path.basename(sys.argv[0])


# print information about username and SVN revision
def getInfo():
    if 'CMSSW_BASE' in os.environ.keys():
        reldir = os.environ['CMSSW_BASE']+"/src"    
        print "--- job_info: begin ---"
        print "SVN revision: " + getSVNrevision(reldir)
    else:
        print "CMSSW enviroment not set"
    print "Username: " + getUsername()
    print "--- job_info: end ---"


# main function
def main():
    args = sys.argv[1:]
    if "-h" in args or "--help" in args:
        usage()
        sys.exit(2)
    getInfo()


if __name__ == "__main__":
   sys.exit(main())
