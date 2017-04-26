#!/bin/bash
#
# Script to update all version numbers in SoPlex
#
# usage: updateversion.sh [<major> <minor> <patch>]"
#

if [[ $# -eq 0 ]]
then
    perl -i -pe 's/#define SOPLEX_APIVERSION        \K(\d+)/$1 + 1/e' src/spxdefines.h
    echo "updated API version:"
    grep SOPLEX_APIVERSION src/spxdefines.h

elif [[ $# -eq 4 ]]
then
    sed -i "s/^VERSION.*/VERSION		:=	$1.$2.$3.$4/" Makefile
    sed -i "s/^VERSION.*/VERSION		=	$1.$2.$3.$4/" Makefile.nmake
    sed -i "s/^ \* @version.*/ \* @version  $1.$2.$3.$4/" doc/xternal.cpp
    sed -i "s/^VERSION=.*/VERSION=\"$1.$2.$3.$4\"/" makedist.sh
    sed -i "s/^#define SOPLEX_VERSION.*/#define SOPLEX_VERSION         $1$2$3/" src/spxdefines.h
    sed -i "s/^#define SOPLEX_SUBVERSION.*/#define SOPLEX_SUBVERSION        $4/" src/spxdefines.h
    echo "new version:"
    grep -e SOPLEX_VERSION -e SOPLEX_SUBVERSION -e SOPLEX_APIVERSION src/spxdefines.h
elif [[ $# -eq 3 ]]
then
    sed -i "s/^VERSION.*/VERSION		:=	$1.$2.$3/" Makefile
    sed -i "s/^VERSION.*/VERSION		=	$1.$2.$3/" Makefile.nmake
    sed -i "s/^ \* @version.*/ \* @version  $1.$2.$3/" doc/xternal.cpp
    sed -i "s/^VERSION=.*/VERSION=\"$1.$2.$3\"/" makedist.sh
    sed -i "s/^#define SOPLEX_VERSION.*/#define SOPLEX_VERSION         $1$2$3/" src/spxdefines.h
    sed -i "s/^#define SOPLEX_SUBVERSION.*/#define SOPLEX_SUBVERSION        0/" src/spxdefines.h
    echo "new version:"
    grep -e SOPLEX_VERSION -e SOPLEX_SUBVERSION -e SOPLEX_APIVERSION src/spxdefines.h
else
    echo "usage: $0 [<major> <minor> <patch>]"
    echo "update all version numbers to the specified one or only increase API version"
    exit 1;
fi
