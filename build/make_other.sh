#!/bin/bash

# stop on error
set -e

# find location
# =============

script="$(readlink -f "$0")"
dir="$(dirname "$script")"


# compile scywalker-report
echo "Building scywalker-report"
cd $dir/../rust/scywalker-report
./install.sh

echo "Finished building genomecomb other tools"
cp $dir/../rust/scywalker-report/target/x86_64-unknown-linux-musl/release/scywalker-report $dir/../extra