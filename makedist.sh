#!/bin/bash -e

VERSION="4.0.0.1"
NAME="soplex-$VERSION"
rm -f $NAME
ln -s . $NAME
rm -f $NAME.tgz

echo "store git hash"
GITHASH=`git describe --always --dirty  | sed 's/^.*-g//'`
echo "#define SPX_GITHASH \"$GITHASH\"" > src/soplex/git_hash.cpp

# Before we create a tarball change the directory and file rights in a command way
echo "adjust file modes"
find ./ -type d -exec chmod 750 {} \;
find ./ -type f -exec chmod 640 {} \;
find ./ -name "*.sh" -exec chmod 750 {} \;
find ./ -name "*.py" -exec chmod 750 {} \;
find ./ -name "*.cmake" -exec chmod 750 {} \;
chmod 750 bin/*

tar -cvzhf $NAME.tgz \
--exclude="*~" \
--exclude=".?*" \
--exclude="*exercise_LP_changes.cpp" \
--exclude="*/local/*" \
--exclude="TODO" \
$NAME/COPYING \
$NAME/INSTALL.md \
$NAME/CHANGELOG \
$NAME/Makefile \
$NAME/check/check.awk \
$NAME/check/check.sh \
$NAME/check/compare.py \
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
$NAME/doc/inc/localfaq.php $NAME/doc/inc/faqtext.txt $NAME/doc/inc/parser.py \
$NAME/doc/builddoc.sh \
$NAME/doc/soplex.dxy \
$NAME/doc/soplexfooter.html \
$NAME/doc/soplexheader.html \
$NAME/doc/soplexlayout.xml \
$NAME/doc/xternal.cpp \
$NAME/make/make* \
$NAME/settings/exact.set \
$NAME/settings/devex.set \
$NAME/settings/steep.set \
$NAME/settings/polish1.set \
$NAME/settings/polish2.set \
$NAME/src/depend* \
$NAME/src/*.h \
$NAME/src/*.cpp \
$NAME/src/soplex/*h \
$NAME/src/soplex/*.cpp \
$NAME/CMakeLists.txt              \
$NAME/settings/default-col.set    \
$NAME/settings/default-row.set    \
$NAME/settings/devex.set          \
$NAME/settings/steep.set          \
$NAME/settings/exact.set          \
$NAME/settings/polish1.set        \
$NAME/settings/polish2.set        \
$NAME/soplex-config.cmake.in      \
$NAME/src/CMakeLists.txt         \
$NAME/check/CMakeLists.txt \
$NAME/cmake/Modules

rm -f $NAME

echo ""
echo "check version numbers in src/soplex/spxdefines.h, doc/xternal.cpp, CMakeLists.txt, Makefile, Makefile.nmake, and makedist.sh ($VERSION):"
grep "VERSION" src/soplex/spxdefines.h
grep "@version" doc/xternal.cpp
grep "SOPLEX_VERSION" CMakeLists.txt
grep "^VERSION" Makefile
grep "^VERSION" Makefile.nmake
grep "^VERSION" makedist.sh
echo "check copyright info in doxygen documentation:"
grep "1996" doc/soplexfooter.html
tail src/soplex/git_hash.cpp
