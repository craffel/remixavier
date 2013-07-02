#!/bin/sh
#
# matlab_arch.sh
#
# Reports machine architecture "as if" it was Matlab's machine(); command

UNAME=`uname -s`
VERS=`uname -r`
MACH=`uname -m`

case $UNAME in

    Darwin)
	case $MACH in
	    x86_64)
		ARCH=MACI64
		arch=maci64
		;;
	    *)
		ARCH=MACI
		arch=maci
		;;
	esac
	;;
    Linux)
	case $MACH in
	    i?86)
		ARCH=GLNX86
		arch=glnx86
		;;
	    x86_64)
		ARCH=GLNXA64
		arch=a64
		;;
	    *)
		ARCH=GLNX_UNK
		arch=glnx
	esac
	;;
    *)
	ARCH=UNKNOWN
	arch=
	;;
esac

if [ $# = 0 ]; then
    echo $ARCH
else
    echo $arch
fi

exit 0
