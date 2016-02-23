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
  #doModel '3B' 'bias'    "$1" 'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
  doModel '4A' 'bias'    "$1" 'random'  'TRUE'  'TRUE'  'FALSE' 'TRUE'
  doModel '5A' 'bias'    "$1" 'average' 'TRUE'  'TRUE'  'TRUE'  'TRUE'
  doModel '6A' 'bias'    "$1" 'random'  'TRUE'  'TRUE'  'TRUE'  'FALSE'
}

doModel '3B' 'bias' "cchf"  'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
doModel '3B' 'bias' "chik"  'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
doModel '3B' 'bias' "deng"  'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
doModel '3B' 'bias' "hat"   'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
doModel '3B' 'bias' "melio" 'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
doModel '3B' 'bias' "nwcl"  'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
doModel '3B' 'bias' "owcl"  'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'
doModel '3B' 'bias' "scrub" 'random'  'TRUE'  'TRUE'  'TRUE'  'TRUE'

doDisease "cchf"
doDisease "chik"
doDisease "deng"
doDisease "hat"
doDisease "melio"
doDisease "nwcl"
doDisease "owcl"
doDisease "scrub"
echo "Done"
date