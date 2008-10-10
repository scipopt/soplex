#!/bin/sh

export PATH="/htcsoft/amd/sles9/gcc/current/bin:$PATH"
export COVDIR="/htcsoft/amd/sles9/coverity/current/bin"
export COVANALYSIS="cov-analysis"
export COVDB="cov-db"
export PRODUCT="soplex"

# checker options: warn about stack size only when
#  - more than 512k are allocated on the stack at once
#  - more than 1M   are allocated on the stack in total
export CHECKEROPTIONS="--checker-option STACK_USE:max_single_base_use_bytes:524288 \
                       --checker-option STACK_USE:max_total_use_bytes:1048576"

#cvs update
$COVDIR/cov-install-gui --password $PRODUCT --datadir $COVDB --product $PRODUCT
make OPT=dbg clean
$COVDIR/cov-build --dir $COVANALYSIS make -s OPT=dbg
$COVDIR/cov-analyze --dir $COVANALYSIS -all $CHECKEROPTIONS
$COVDIR/cov-commit-defects --datadir $COVDB --product $PRODUCT --user admin --dir $COVANALYSIS
$COVDIR/cov-start-gui --datadir $COVDB --port 9402

