#!/usr/bin/env bash

# very simple evaluation script

OUTFILE=$1

optimal=$(expr $(grep -c "SoPlex status       : problem is solved \[optimal\]" $OUTFILE) / 2)
total=$(grep -c "=ready=" $OUTFILE)

echo
echo 'Summary:'
echo
echo 'total  : '$total
echo 'optimal: '$optimal
echo
echo 'please install python for better statistics and evaluation'
