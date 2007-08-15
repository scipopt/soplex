#!/bin/sh
tar -cvzhf soplex-1.3.1.tgz --files-from=soplex-1.3.1/distrib.list --exclude="*CVS*" --exclude="*~" --exclude=".?*" --exclude="*exercise_LP_changes.cpp" --exclude="*/local/*" 
