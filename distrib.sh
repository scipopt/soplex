#!/bin/sh
tar -cvzhf soplex-1.2.2-beta2.tar.gz --files-from=soplex-1.2.2/distrib.list --exclude="*CVS*" --exclude="*~" --exclude=".?*"
