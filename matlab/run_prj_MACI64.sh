#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
# Modified by dpwe@ee.columbia.edu to assume default MCR install location
#   run_renoiser_prj_MACI64.sh - for Matlab R2010b on 64-bit Mac Intel
# $Header: /Users/dpwe/docs/grants/2010-01-DARPA-RATS/code/snreval/RCS/run_snreval.sh,v 1.1 2010/12/14 16:32:10 dpwe Exp dpwe $

exe_name=$0
exe_dir=`dirname "$0"`

MCRROOT=/Applications/MATLAB_R2010b.app

MWE_ARCH="maci64" ;
DYLD_LIBRARY_PATH=.:${MCRROOT}/runtime/${MWE_ARCH} ;
DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/bin/${MWE_ARCH} ;
DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/sys/os/${MWE_ARCH} ;
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export DYLD_LIBRARY_PATH;
export XAPPLRESDIR;
#echo DYLD_LIBRARY_PATH is ${DYLD_LIBRARY_PATH};

args=
while [ $# -gt 0 ]; do
    token=`echo "$1" | sed -e "s/[][\' ()&;^]/\\\\\&/g"`   # Add blackslash before each blank and backslash
    args="${args} ${token}" 
    shift
done
#prjname=`echo $0 | sed -e "s@.*/@@" -e 's@\.[^.]*$@@' -e "s/^run_//" -e "s/_maci64//"`
prjname=`echo $0 | sed -e "s@.*/@@" -e 's@\.[^.]*$@@' -e "s/^run_//"`
#echo "args=$args"
eval "${exe_dir}"/${prjname}_${MWE_ARCH}.app/Contents/MacOS/${prjname}_${MWE_ARCH} $args

exit

