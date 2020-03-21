#!/bin/bash

set -xeuo pipefail

ARGS="$@"
# directory handed in over symlink from corresponding fetch step
SRC=fetch/mini-nest
for file in $SRC/main.cpp $SRC/stopwatch.cpp $SRC/simulate.cpp; do \
  $compile $ARGS -c $SRC/main.cpp;
done
$link main.o stopwatch.o simulate.o -o main
