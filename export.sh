#!/usr/bin/env bash
# 
#  bump-version.sh
#  ObjectModel
#  
#  Created by Alexander Rudy on 2011-10-17.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 


EXE=$0
BRANCH="master"
EXCLUDE="exclusions.txt"
USAGE="
Usage for bump-version

	$EXE path/to/destination/directory
	
	Will export branch $BRANCH to the destination directory, excluding files listed in $EXCLUDE
"

if [ "$#" -ne 1 ]
then
	echo "$USAGE"
	exit
fi

DESTINATION=$1

echo "Exporting branch $BRACNH to $DESTINATION"

git archive $BRANCH | tar x -C $DESTINATION

for file in `cat $EXCLUDE`
do
	rm "$DESTINATION/$file"
	echo "Removing $DESTINATION/$file"
done

echo "Done."