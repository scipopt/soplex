#!/bin/sh
tar -cvzhf soplex-1.3.0.tar.gz --files-from=soplex-1.3.0/distrib.list --exclude="*CVS*" --exclude="*~" --exclude=".?*" --exclude="*exercise_LP_changes.cpp" --exclude="*/local/*" 
