#!/bin/bash
#
# This bash script updates all copyrights in the SoPlex files and posted
# those files which contain a copyright which do not have the right format
#
# You just have to run this script. There is nothing to adjust.
# The correct year is detected through the 'date' function
#
# Note that not all files (usually scripts) contain a copyright. A copyright is only
# needed for those files which are part of a SoPlex distribution (see scripts/makedist.sh)
#
# USAGE: ./scripts/updatedates.sh

NEWYEAR=`date +"%Y"`
LASTYEAR=`expr $NEWYEAR - 1`

echo ""
echo "NEWYEAR  = $NEWYEAR"
echo "LASTYEAR = $LASTYEAR"

DIRECTORIES=(check doc src src/soplex tests extra)
EXTENSIONS=(sh awk h c hpp cpp html)
EXTRAFILES=(Makefile make/make.install make/make.detecthost Makefile.nmake)

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
      COUNT1=`grep -c "1996-$NEWYEAR" $FILE`
      COUNT2=`grep -c "1996-$LASTYEAR" $FILE`

      if test "$COUNT2" != 0
      then
          if test "$COUNT1" == "$COUNT2"
          then
              echo "- $FILE PARTIALLY UP TO DATE, updating completely"
          else
              echo "- $FILE updated"
          fi

          sed -i 's!1996-'$LASTYEAR'!1996-'$NEWYEAR'!g' $FILE
      else
          echo "- $FILE already up to date"
          continue
      fi
            fi
        done
    done
done

echo ""
echo "Updating additional files.."

for FILE in ${EXTRAFILES[@]}
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
      COUNT1=`grep -c "1996-$NEWYEAR" $FILE`
      COUNT2=`grep -c "1996-$LASTYEAR" $FILE`

      if test "$COUNT2" != 0
      then
          if test "$COUNT1" == "$COUNT2"
          then
              echo "- $FILE PARTIALLY UP TO DATE, updating completely"
          else
              echo "- $FILE updated"
          fi

          sed -i 's!1996-'$LASTYEAR'!1996-'$NEWYEAR'!g' $FILE
      else
          echo "- $FILE already up to date"
          continue
      fi
    fi
done

sed -i 's!2002-'$LASTYEAR'!2002-'$NEWYEAR'!g' LICENSE
