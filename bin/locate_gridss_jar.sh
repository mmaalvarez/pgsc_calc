#!/bin/bash

gridssShellPath=`which gridss`
set -o pipefail

SOURCE="$gridssShellPath"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# FIX: Make output deterministic
find $DIR -maxdepth 1 -name "gridss*.jar" -o -name "picard*.jar" | sort | head -n 1
