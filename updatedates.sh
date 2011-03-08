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
# $Id: updatedates.sh,v 1.6 2011/02/24 13:53:37 bzfgleix Exp $

NEWYEAR=`date +"%Y"`
LASTYEAR=`expr $NEWYEAR - 1`

echo ""
echo "NEWYEAR  = $NEWYEAR"
echo "LASTYEAR = $LASTYEAR"

DIRECTORIES=(check doc src tests extra)
EXTENSIONS=(sh awk h c hpp cpp)
EXTRAFILES=(Makefile INSTALL extra/lpconv.cpp extra/lpstat.cpp)

for DIRECTORY in ${DIRECTORIES[@]}
do
  echo ""
  echo "Updating directory $DIRECTORY.."

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
	  COUNT1=`grep -c "1996\-$NEWYEAR" $FILE`
	  COUNT2=`grep -c "1996\-$LASTYEAR" $FILE`

	  if test "$COUNT2" != 0
	  then
	      if test "$COUNT1" == "$COUNT2"
	      then
		  echo "- $FILE PARTIALLY UP TO DATE, updating completely"
	      else
		  echo "- $FILE updated"
	      fi
  
	      mv $FILE $FILE.olddate
	      sed 's!1996-'$LASTYEAR'!1996-'$NEWYEAR'!g' $FILE.olddate > $FILE
	      rm $FILE.olddate
	  else
	      echo "- $FILE already up to date"
	      continue
	  fi
      fi
    done
  done
done

