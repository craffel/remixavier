#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
# Modified by dpwe@ee.columbia.edu to assume default MCR install location
#  run_renoiser_prj_GLNXA64.sh - for Matlab R2007 (v77) on 64-bit x86_64 Linux
# $Header: /Users/dpwe/docs/grants/2010-01-DARPA-RATS/code/snreval/RCS/run_snreval.sh,v 1.1 2010/12/14 16:32:10 dpwe Exp dpwe $

exe_name=$0
exe_dir=`dirname "$0"`

MCRROOT=/opt/MATLAB/MATLAB_Compiler_Runtime/v714

MWE_ARCH="glnxa64" ;
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/${MWE_ARCH} ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/${MWE_ARCH} ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/${MWE_ARCH};
        MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export LD_LIBRARY_PATH;
export XAPPLRESDIR;

args=
while [ $# -gt 0 ]; do
    token=`echo "$1" | sed -e "s/[][\' ()&;^"\$"]/\\\\\&/g"`   # Add blackslash before each blank and backslash
    args="${args} ${token}" 
    shift
done
prjname=`echo $0 | sed -e "s@.*/@@" -e 's@\.[^.]*$@@' -e "s/^run_//"`
eval "${exe_dir}"/${prjname}_prj $args

exit

