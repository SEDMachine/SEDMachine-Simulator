#!/usr/bin/env bash
# 
#  get-AstroObject.sh
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2011-11-23.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

CMD=$0
ARG=$1
URL="https://github.com/alexrudy/AstroObject/tarball/master"
GET="curl -L"
TAR="tar"
EXE=${CMD##*/}
LOC=${CMD%%*/}
CWD=`pwd`

USAGE="
$EXE - Fetch the latest release of AstroObject

Usage:
	$EXE
	$EXE -h

$EXE uses $RSYNC to copy from a target to a backup destination.
	The version of $RSYNC used is $RSYNCV

The local directory is currently set to:
LOCAL: $LOCAL

The targets are located in the remote directory:
SERVER: $SERVER
DIRECTORY: $REMOTE

The targets are

	all			runs all sync requirements
	images		final images
	partials	debugging images and files
	logs		system logs
	caches		cache contents
"

if [ $# -gt 0 ]
then
	echo -e "$USAGE"
	exit
fi
if [ $CWD != $LOC ]
then
	echo -e "Must be in identical Directory"
	exit
fi

if [ -d "$CWD/AstroObject/.git" ]
then
	echo -e "Will not overwrite a live git repo!"
	exit
fi

$GET $URL | $TAR xvz -C "$CWD"


rm -r "$CWD/AstroObject"

mkdir "$CWD/AstroObject"

DEST=`ls -1 $CWD | grep alexrudy-AstroObject`

echo "Moving from $DEST to AstroObject"

mv "$CWD/$DEST/"* "$CWD/AstroObject"


rm -r "$CWD/$DEST"








