#!/usr/bin/env bash
set -e
function doModel {
  if [ -d "temp" ]; then
    rm -rf "temp"
  fi
  echo "$3_$1"
  echo -e "source('runTest.R')\n" > go.R
  echo -e "runTest(name='$3_$1', mode='$2', disease='$3', admin_extract_mode='$4', crop_bias=$5, filter_bias=$6, use_weights=$7, use_temporal_covariates=$8)\n" >> go.R
  date
  R --no-save --slave -f go.R &> "$3_$1.log"
  rm go.sh
  rm -rf "temp"
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
doDisease "owcl"  "FALSE"
doDisease "scrub" "TRUE"
echo "Done"
date