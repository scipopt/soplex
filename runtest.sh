#!/bin/sh
ARCH=`uname -m | sed -e s/sun../sparc/ -e s/i.86/x86/ -e s/IP../mips/ -e s/9000..../hppa/`
OSTYPE=`uname -s | tr A-Z a-z`
case $OSTYPE in
linux)
   case $ARCH in
   x86)
      gmake COMP=gnu    OPT=opt clean
      gmake COMP=gnu    OPT=opt
      gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
      ;;
   alpha)
      gmake COMP=compaq OPT=opt clean 
      gmake COMP=compaq OPT=opt
      gmake COMP=compaq OPT=opt check >/dev/null 2>&1 &
      gmake COMP=gnu    OPT=opt clean
      gmake COMP=gnu    OPT=opt
      gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
      ;;
   esac
osf1)
   gmake COMP=compaq OPT=opt clean 
   gmake COMP=compaq OPT=opt
   gmake COMP=compaq OPT=opt check >/dev/null 2>&1 &
   gmake COMP=gnu    OPT=opt clean
   gmake COMP=gnu    OPT=opt
   gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
   ;;
sunos)
   gmake COMP=sun    OPT=opt clean 
   gmake COMP=sun    OPT=opt
   gmake COMP=sun    OPT=opt check >/dev/null 2>&1 &
#   gmake COMP=gnu    OPT=opt clean
#   gmake COMP=gnu    OPT=opt
#   gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
   ;;
hp-ux)
   gmake COMP=sun    OPT=opt clean 
   gmake COMP=sun    OPT=opt
   gmake COMP=sun    OPT=opt check >/dev/null 2>&1 &
   ;;
irix)
   gmake COMP=sgi    OPT=std clean 
   gmake COMP=sgi    OPT=std
   gmake COMP=sgi    OPT=std check >/dev/null 2>&1 &
#   gmake COMP=gnu    OPT=opt clean
#   gmake COMP=gnu    OPT=opt
#   gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
   ;;
esac



