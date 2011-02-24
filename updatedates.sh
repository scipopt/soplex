#!/bin/bash
#
# This bash script updates all copyrights in the SoPlex files and posted
# those files which contain a copyright which do not have the right format
#
# You just have to run this script. There is nothing to adjust. 
# The correct year is detected through the 'date' function 
#
# Note that not all files (usually scripts) contain a copyright. A copyright is only 
# needed for those files which are part of a SoPlex distribution (see makedist.sh)
#
# $Id: updatedates.sh,v 1.5 2011/02/24 11:52:52 bzfgleix Exp $

NEWYEAR=`date +"%Y"`
LASTYEAR=`expr $NEWYEAR - 1`

DIRECTORIES=(check doc src tests extra)
EXTENSIONS=(sh awk h c hpp cpp)
EXTRAFILES=(Makefile INSTALL extra/lpstat.cpp extra/lpconv.cpp)

for DIRECTORY in ${DIRECTORIES[@]}
do
  for EXTENSION in ${EXTENSIONS[@]}
  do
    for FILE in $DIRECTORY/*.$EXTENSION
    do
      if test -f $FILE
      then
	  # check if the file has a correct old copyright 
	  COUNT=`grep -c 1997- $FILE`

	  if test "$COUNT" != 0
	  then
	      # post those files which have a wrong old copyright
	      echo "COPYRIGHT ERROR --------------------> $FILE"
	      grep "1997-" $FILE
	      continue
	  fi

	  # check if the file has already the new copyright 
	  COUNT=`grep -c -$NEWYEAR $FILE`

	  if test "$COUNT" != 0
	  then
	      continue
	  else
	      echo $FILE
  
	      mv $FILE $FILE.olddate
	      sed 's!1996-'$LASTYEAR'!1996-'$NEWYEAR'!g' $FILE.olddate > $FILE
	      rm $FILE.olddate
	  fi
      fi
    done
  done
done

