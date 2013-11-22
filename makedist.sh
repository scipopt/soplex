#!/bin/sh

VERSION="1.7.0.5"
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
$NAME/check/instances/* \
$NAME/doc/soplex.dxy \
$NAME/doc/soplex.css \
$NAME/doc/soplexfooter.html \
$NAME/doc/soplexheader.html \
$NAME/doc/soplexlayout.xml \
$NAME/doc/xternal.cpp \
$NAME/doc/inc/faq.inc \
$NAME/doc/inc/faqcss.inc \
$NAME/make/make* \
$NAME/src/depend* \
$NAME/src/*h \
$NAME/src/*cpp
rm -f $NAME

echo ""
echo "check version numbers in src/spxdefines.h, doc/xternal.cpp, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/spxdefines.h
grep "@version" doc/xternal.cpp
grep "^VERSION" Makefile
grep "^VERSION" makedist.sh
echo "check copyright info in doxygen documentation:"
grep "2003" doc/soplexfooter.html
tail src/git_hash.cpp
