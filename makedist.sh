#!/bin/sh

VERSION="2.0.1.2"
NAME="soplex-$VERSION"
rm -f $NAME
ln -s . $NAME
rm -f $NAME.tgz

# compile to create the correct GiTHash
make githash

# Before we create a tarball change the directory and file rights in a command way
echo adjust file modes
find ./ -type d -exec chmod 750 {} \;
find ./ -type f -exec chmod 640 {} \;
find ./ -name "*.sh" -exec chmod 750 {} \;
find ./ -name "*.py" -exec chmod 750 {} \;
chmod 750 bin/*

tar -cvzhf $NAME.tgz \
--exclude="*CVS*" \
--exclude="*~" \
--exclude=".?*" \
--exclude="*exercise_LP_changes.cpp" \
--exclude="*/local/*" \
--exclude="TODO" \
$NAME/COPYING \
$NAME/INSTALL \
$NAME/CHANGELOG \
$NAME/Makefile \
$NAME/check/check.sh \
$NAME/check/check_legacy.sh \
$NAME/check/evaluation.py \
$NAME/check/evaluation.sh \
$NAME/check/instances/* \
$NAME/check/test.sh \
$NAME/check/testset/netlib.test \
$NAME/check/testset/quick.test \
$NAME/check/testset/mittelmann.test \
$NAME/check/testset/infeas.test \
$NAME/check/testset/netlib.solu \
$NAME/check/testset/quick.solu \
$NAME/check/testset/mittelmann.solu \
$NAME/check/testset/infeas.solu \
$NAME/doc/soplex.dxy \
$NAME/doc/soplex.css \
$NAME/doc/soplexfooter.html \
$NAME/doc/soplexheader.html \
$NAME/doc/soplexlayout.xml \
$NAME/doc/xternal.cpp \
$NAME/doc/inc/faq.inc \
$NAME/doc/inc/faqcss.inc \
$NAME/make/make* \
$NAME/settings/default.set \
$NAME/settings/devex.set \
$NAME/settings/steep.set \
$NAME/src/depend* \
$NAME/src/*h \
$NAME/src/*cpp
rm -f $NAME

echo ""
echo "check version numbers in src/spxdefines.h, doc/xternal.cpp, Makefile, Makefile.nmake, and makedist.sh ($VERSION):"
grep "VERSION" src/spxdefines.h
grep "@version" doc/xternal.cpp
grep "^VERSION" Makefile
grep "^VERSION" Makefile.nmake
grep "^VERSION" makedist.sh
echo "check copyright info in doxygen documentation:"
grep "2003" doc/soplexfooter.html
tail src/git_hash.cpp
