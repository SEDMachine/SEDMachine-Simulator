#!/usr/bin/env bash
# 
#  get-images.sh
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2011-11-23.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 
CMD=$0
DEST=$1
REMOTE="~/SEDM"
SERVER="sedmachine"
LOCAL="$HOME/Dropbox/SEDM/SEDM-Dev"
RSYNC='rsync'
RSYNCV=`rsync --version | head -n1`
EXE=${CMD##*/}

USAGE="
$EXE - Collection of files on $SERVER using $RSYNC

Usage:
	$EXE [command]
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

case $DEST in
    all )
        NAME="all"
		echo "all not implemented -- sorry!"
		exit
        ;;
    images )
        NAME="images"
        MASTER="$REMOTE/Images/"
        DESTINATION="$LOCAL/Images/"
        ;;
	partials )
		NAME="partials"
		MASTER="$REMOTE/Partials/"
		DESTINATION="$LOCAL/Partials/"
		;;
	logs )
		NAME="logs"
		MASTER="$REMOTE/Logs/"
		DESTINATION="$LOCAL/Logs/"
		;;
	caches )
		NAME="caches"
		MASTER="$REMOTE/Caches/"
		DESTINATION="$LOCAL/Caches/"
		;;
	*)
		echo -e "$USAGE"
		exit
		;;
esac

mkdir -p "$DESTINATION"
echo "Backing up $NAME"
echo "DRY RUN from $SERVER:$MASTER to $DESTINATION"
$RSYNC -avn --progress "$SERVER:$MASTER" "$DESTINATION"
echo "Continue with $NAME RSYNC from $SERVER:$MASTER to $DESTINATION? (y/n)"
read cont
if [ $cont != "y" ]
then
    echo "Cancelling backup from $SERVER:$MASTER to $DESTINATION"
	exit
fi
$RSYNC -av "$SERVER:$MASTER" "$DESTINATION"
echo ""
echo "Backed Up from $SERVER:$MASTER to $DESTINATION"
echo "Thanks for backing up your $NAME!"
