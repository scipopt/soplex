#!/bin/sh
for i in Makefile INSTALL check/*.sh check/*.awk doc/*.cpp src/*.cpp src/*.h tests/*.cpp extra/*.cpp
do
if [ -f $i ]
then
echo $i
mv $i $i.olddate
sed 's!1997-2008!1997-2009!g' $i.olddate > $i
fi
done
