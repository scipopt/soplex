#!/bin/sh
VERSION="1.4.1"
NAME="soplex-$VERSION"
rm -f $NAME
ln -s . $NAME
rm -f $NAME.tgz
tar -cvzhf $NAME.tgz \
--exclude="*CVS*" \
--exclude="*~" \
--exclude=".?*" \
--exclude="*exercise_LP_changes.cpp" \
--exclude="*/local/*" \
$NAME/COPYING \
$NAME/README \
$NAME/CHANGELOG \
$NAME/Makefile \
$NAME/check/netlib.test \
$NAME/check/quick.test \
$NAME/check/mittelmann.test \
$NAME/check/infeas.test \
$NAME/check/netlib.solu \
$NAME/check/quick.solu \
$NAME/check/mittelmann.solu \
$NAME/check/infeas.solu \
$NAME/check/check.sh \
$NAME/check/check.awk \
$NAME/doc/soplex.dxy \
$NAME/doc/xternal.cpp \
$NAME/doc/inc/faq.inc \
$NAME/doc/inc/faqcss.inc \
$NAME/make/* \
$NAME/src/*
rm -f $NAME
echo ""
echo "check version numbers in src/spxdefines.h, doc/xternal.cpp, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/spxdefines.h
grep "@version" doc/xternal.cpp
grep "^VERSION" Makefile
