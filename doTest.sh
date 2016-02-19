#!/usr/bin/env bash
set -e
function doModel {
	if [ -d "temp" ]; then
		rm -rf "temp"
	fi
	echo "$3_$1"
	date
	R --no-save --slave -e "source('runTest.R'); runTest('$3_$1', '$2', '$3', '$4', '$5', '$6', '$7', '$8');" &> "$3_$1.log"
	if [ -d "temp" ]; then
		rm -rf "temp"
	fi
}

function doDisease {
	doModel '1A' 'bhatt'   "$1" 'random'  'FALSE' 'FALSE' 'TRUE'  'TRUE'
	doModel '1B' 'bias'    "$1" 'random'  'FALSE' 'FALSE' 'TRUE'  'TRUE'
	doModel '1C' 'uniform' "$1" 'random'  'FALSE' 'FALSE' 'TRUE'  'TRUE'
	doModel '2C' 'bias'    "$1" 'random'  'TRUE'  'FALSE' 'TRUE'  'TRUE'
	doModel '2D' 'uniform' "$1" 'random'  'TRUE'  'FALSE' 'TRUE'  'TRUE'
	doModel '3B' 'bias'    "$1" 'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
	doModel '4A' 'bias'    "$1" 'random'  'TRUE'  'TRUE'  'FALSE' 'TRUE'
	doModel '5A' 'bias'    "$1" 'average' 'TRUE'  'TRUE'  'TRUE'  'TRUE'
	if [ "$2" == "TRUE" ]; then
		doModel '6A' 'bias' '$1' 'random' 'TRUE'  'TRUE'  'TRUE'  'FALSE'
	fi
}

doDisease "cchf"  "TRUE"
doDisease "chik"  "FALSE"
doDisease "deng"  "FALSE"
doDisease "hat"   "TRUE"
doDisease "melio" "TRUE"
doDisease "nwcl"  "FALSE"
doDisease "nwvl"  "FALSE"
doDisease "owcl"  "FALSE"
doDisease "owvl"  "FALSE"
doDisease "scrub" "TRUE"
echo "Done"
date