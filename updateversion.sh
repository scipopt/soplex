#!/bin/bash
#
# Script to update all version numbers in SoPlex
#

if [ $# -eq 4 ]
then
    sub=".$4"
    subzero="        $4"
elif [ $# -eq 3 ]
then
    sub=""
    subzero="        0"
else
    echo "usage: $0 <major_1> <major_2> <major_3> [<minor>]"
    echo "update all version numbers to the specified one"
    echo "example dev:     $0 3 0 1 2"
    echo "example release: $0 3 0 2"
    exit 1;
fi

sed -i "s/^VERSION.*/VERSION		:=	$1.$2.$3$sub/" Makefile
sed -i "s/^VERSION.*/VERSION		=	$1.$2.$3$sub/" Makefile.nmake
sed -i "s/^ \* @version.*/ \* @version  $1.$2.$3$sub/" doc/xternal.cpp
sed -i "s/^VERSION=.*/VERSION=\"$1.$2.$3$sub\"/" makedist.sh
sed -i "s/^#define SOPLEX_VERSION.*/#define SOPLEX_VERSION         $1$2$3/" src/spxdefines.h
sed -i "s/^#define SOPLEX_SUBVERSION.*/#define SOPLEX_SUBVERSION$subzero/" src/spxdefines.h
