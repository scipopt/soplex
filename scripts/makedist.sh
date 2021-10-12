#!/bin/bash -e

# create tarball for release
# usage: ./scripts/makedist.sh

VERSION="5.0.2.4"
NAME="soplex-$VERSION"
rm -f $NAME.tgz
rm -f $NAME.tar

echo "store git hash"
GITHASH=`git describe --always --dirty  | sed 's/^.*-g//'`
echo "#define SPX_GITHASH \"$GITHASH\"" > src/soplex/git_hash.cpp

# Before we create a tarball change the directory and file rights in a command way
echo "adjust file modes"
git ls-files | xargs dirname | sort -u | xargs chmod 750
git ls-files | xargs chmod 640
git ls-files "*.sh" "*.py" | xargs chmod 750

# pack files tracked by git and append $NAME to the front
git ls-files -c | xargs tar --transform "s|^|${NAME}/|" -cvhf $NAME.tar \
    --exclude="*~" \
    --exclude=".*" \
    --exclude="check/*cluster*" \
    --exclude="check/make_solu*" \
    --exclude="check/testset/*" \
    --exclude="extra/*" \
    --exclude="lcov/*" \
    --exclude="license/*" \
    --exclude="lint/*" \
    --exclude="make/local/*" \
    --exclude="paper/*" \
    --exclude="scripts/makedist.sh" \
    --exclude="scripts/update*.sh" \
    --exclude="tests/*.cpp" \
    --exclude="tests/Makefile" \
    --exclude="web/*"

# append additional files that were excluded before
tar --transform "s|^|${NAME}/|" -rvf $NAME.tar \
    check/testset/quick.test \
    check/testset/quick.solu \
    src/soplex/git_hash.cpp

# compress the archive
gzip -c $NAME.tar > $NAME.tgz

# remove temporary archive
rm -f $NAME.tar

echo ""
echo "check version numbers ($VERSION):"
grep -H "VERSION" src/soplex/spxdefines.h
grep -H "@version" doc/xternal.cpp
grep -H "SOPLEX_VERSION" CMakeLists.txt
grep -H "^VERSION" Makefile
grep -H "^VERSION" scripts/makedist.sh
echo "check copyright info in doxygen documentation:"
grep "1996" doc/soplexfooter.html
tail src/soplex/git_hash.cpp
