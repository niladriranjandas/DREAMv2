#!/bin/bash

proteinname="$1"

./convertToXplor.sh "$proteinname"

./driveCode.sh "$proteinname"

./callBadGap.sh "$proteinname"
