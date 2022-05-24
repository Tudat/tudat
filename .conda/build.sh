#!/usr/bin/env bash

mkdir build

cd build

cmake \
    -DCMAKE_CXX_STANDARD=14 \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_PREFIX_PATH=$PREFIX \
    -DPREFIX=$PREFIX \
    -DTUDAT_CONDA_BUILD=on \
    -DTUDAT_BUILD_STATIC_LIBRARY=on \
    -DTUDAT_BUILD_TUDAT_TUTORIALS=off \
    -DTUDAT_BUILD_WITH_SOFA_INTERFACE=on \
    -DTUDAT_BUILD_WITH_SPICE_INTERFACE=on \
    -DTUDAT_INSTALL=on \
    -DTUDAT_TEST_INSTALL=off \
    ..

make -j2

ctest

make install
