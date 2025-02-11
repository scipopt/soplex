#!/bin/bash
#
# Script to update all version numbers in SoPlex
#
# USAGE: ./scripts/updateversion.sh [<major> <minor> <patch>]"
#

if [[ $1 == "-a" ]]
then
    perl -i -pe 's/#define SOPLEX_APIVERSION        \K(\d+)/$1 + 1/e' src/soplex/spxdefines.h
    echo "updated API version:"
    grep SOPLEX_APIVERSION src/soplex/spxdefines.h

elif [[ $# -eq 3 ]]
then
    sed -i "s/^VERSION.*/VERSION		:=	$1.$2.$3/" Makefile
    sed -i "s/^ \* @version.*/ \* @version  $1.$2.$3/" doc/xternal.cpp
    sed -i "s/^VERSION=.*/VERSION=\"$1.$2.$3\"/" scripts/makedist.sh
    sed -i "s/^#define SOPLEX_VERSION.*/#define SOPLEX_VERSION         $1$2$3/" src/soplex/spxdefines.h
    sed -i "s/set(SOPLEX_VERSION_MAJOR.*/set(SOPLEX_VERSION_MAJOR $1)/" CMakeLists.txt
    sed -i "s/set(SOPLEX_VERSION_MINOR.*/set(SOPLEX_VERSION_MINOR $2)/" CMakeLists.txt
    sed -i "s/set(SOPLEX_VERSION_PATCH.*/set(SOPLEX_VERSION_PATCH $3)/" CMakeLists.txt

    echo "new version:"
    grep -e SOPLEX_VERSION -e SOPLEX_APIVERSION src/soplex/spxdefines.h
    grep "set(SOPLEX_VERSION" CMakeLists.txt
else
    echo "usage:"
    echo ""
    echo "$0 <major> <minor> <patch> [<sub>]"
    echo " -- update all version numbers to the specified one"
    echo ""
    echo "$0 -a"
    echo " -- only increase API version"
    echo ""
    echo "current version (from spxdefines.h): "
    grep -e SOPLEX_VERSION -e SOPLEX_APIVERSION src/soplex/spxdefines.h
    exit 1;
fi
