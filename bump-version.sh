#!/usr/bin/env bash
# 
#  bump-version.sh
#  ObjectModel
#  
#  Created by Alexander Rudy on 2011-10-17.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 


EXE=$0
USAGE="
Usage for bump-version

	$EXE x.x.x
	
	change the version number to x.x.x

"

if [ "$#" -ne 1 ]
then
	echo "$USAGE"
	exit
fi

SELECTREGEX="[0-9a-zA-Z\.\+]+"
DIR="SEDMachine"
VSPECFILE="VERSION"
AVSPECFILE="$DIR/VERSION"

VERSION=`cat $VSPECFILE`

echo "Version is currently $VERSION, changing to $1"

echo "$1" > $VSPECFILE
echo "$1" > $AVSPECFILE

VERSION=`cat $VSPECFILE`

echo "New Version $VERSION"

echo "Manipulating Python (.py) files"
files=`find . -name '*.py'`

for file in $files
do
	sed -i '' -Ee "s/# +Version $SELECTREGEX/#  Version $VERSION/" $file
	sed -i '' -Ee "s/__version__ += +\'$SELECTREGEX\'/__version__ = \'$VERSION\'/" $file
	sed -i '' -Ee "s/version\s+=\s+[\'\"]$SELECTREGEX[\'\"]/version = \'$VERSION\'/" $file
	sed -i '' -Ee "s/release += +\'$SELECTREGEX\'/version = \'$VERSION\'/" $file
	
	echo "  Changed Version to $VERSION in $file"
done

files=`find . -name '*.md'`

echo "Manipulating Markdown (.md) files"
for file in $files
do
	sed -i '' -Ee "s/ +Version $SELECTREGEX/  Version $VERSION/" $file
	echo "  Changed Version to $VERSION in $file"
done

echo "Done."
