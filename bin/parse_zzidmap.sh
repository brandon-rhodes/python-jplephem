#!/bin/bash

awk '/^ *BLTCOD/ {n=$NF} /^ *BLTNAM/ {sub(".*'\''", "");print}' "$@"
