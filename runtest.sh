#!/bin/sh
# $Id: runtest.sh,v 1.5 2002/01/12 22:02:44 bzfkocht Exp $
ARCH=`uname -m | sed -e s/sun../sparc/ -e s/i.86/x86/ -e s/IP../mips/ -e s/9000..../hppa/`
OSTYPE=`uname -s | tr A-Z a-z`
case $OSTYPE in
linux)
   case $ARCH in
   x86)
      cd /optimi/kombadon/bzfkocht/soplex
      gmake COMP=gnu    OPT=opt clean
      gmake COMP=gnu    OPT=opt
      gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
      ;;
   alpha)
      cd /optimi/kombadon/bzfkocht/soplex
      gmake COMP=compaq OPT=opt clean 
      gmake COMP=compaq OPT=opt
      gmake COMP=compaq OPT=opt check >/dev/null 2>&1 &
      gmake COMP=gnu    OPT=opt clean
      gmake COMP=gnu    OPT=opt
      gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
      ;;
   esac
   ;;
osf1)
   PATH=$PATH:/client/bin
   cd /optimi/kombadon/bzfkocht/soplex
   gmake COMP=compaq OPT=opt clean 
   gmake COMP=compaq OPT=opt
   gmake COMP=compaq OPT=opt check >/dev/null 2>&1 &
   gmake COMP=gnu    OPT=opt clean
   gmake COMP=gnu    OPT=opt
   gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
   ;;
sunos)
   cd /optimi/kombadon/bzfkocht/soplex
   gmake COMP=sun    OPT=opt clean 
   gmake COMP=sun    OPT=opt
   gmake COMP=sun    OPT=opt check >/dev/null 2>&1 &
#   gmake COMP=gnu    OPT=opt clean
#   gmake COMP=gnu    OPT=opt
#   gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
   ;;
hp-ux)   
   PATH=$PATH:/usr/local/bin:/client/bin
   CVS_RSH=ssh
   export CVS_RSH
   cd $HOME/soplex
   cvs update
   gmake COMP=hp     OPT=opt clean 
   gmake COMP=hp     OPT=opt
   gmake COMP=hp     OPT=opt check >/dev/null 2>&1 &
   ;;
irix)
   PATH=$PATH:/usr/local/bin:/client/bin
   CVS_RSH=ssh
   export CVS_RSH
   cd $HOME/soplex
   cvs update
   gmake COMP=sgi    OPT=std clean 
   gmake COMP=sgi    OPT=std
   gmake COMP=sgi    OPT=std check >/dev/null 2>&1 &
#   gmake COMP=gnu    OPT=opt clean
#   gmake COMP=gnu    OPT=opt
#   gmake COMP=gnu    OPT=opt check >/dev/null 2>&1 &
   ;;
esac



