#!/usr/bin/env bash
this=$(dirname "${BASH_SOURCE[0]}")
mkdir -p "$this"/build
pushd "$this"/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD" \
  -DENABLE_GPU=ON ..
make -j install
popd
