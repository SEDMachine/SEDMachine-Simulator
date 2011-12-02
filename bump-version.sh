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

VSPECFILE="VERSION"

VERSION=`cat $VSPECFILE`

echo "Version is currently $VERSION, changing to $1"

echo "$1" > $VSPECFILE

VERSION=`cat $VSPECFILE`

echo "New Version $VERSION"

echo "Manipulating Python (.py) files"
files=`find *.py`

for file in $files
do
	sed -i '' -Ee "s/# +Version [0-9\.]+/#  Version $VERSION/" $file
	echo "  Changed Version to $VERSION in $file"
done

files=`find *.md`

echo "Manipulating Markdown (.md) files"
for file in $files
do
	sed -i '' -Ee "s/ +Version [0-9\.]+/  Version $VERSION/" $file
	echo "  Changed Version to $VERSION in $file"
done

echo "Manipulating Special Files:"
sed -i '' -Ee "s/__version__ += +\'[0-9\.]+\'/__version__ = \'$VERSION\'/" '__init__.py'
echo "  Changed __init__.py version variable to $VERSION"
sed -i '' -Ee "s/version += +\'[0-9\.]+\'/version = \'$VERSION\'/" 'Docs/conf.py'
echo "  Changed Sphinyx conf.py version variable to $VERSION"
sed -i '' -Ee "s/release += +\'[0-9\.]+\'/release = \'$VERSION\'/" 'Docs/conf.py'
echo "  Changed Sphinyx conf.py release variable to $VERSION"

echo "Done."