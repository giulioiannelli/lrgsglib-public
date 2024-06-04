#!/bin/sh

# Determine the project root directory based on the location of this script
LRGSG_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Export directories
export LRGSG_ROOT
export LRGSG_ROOT
export LRGSG_BUILD="${LRGSG_ROOT}/build"
export LRGSG_DATA="${LRGSG_ROOT}/data"
export LRGSG_IPYNB="${LRGSG_ROOT}/ipynb"
export LRGSG_SRC="${LRGSG_ROOT}/src"
export LRGSG_LIB="${LRGSG_SRC}/lrgsglib"
export LRGSG_LIB_CCORE="${LRGSG_LIB}/Ccore"
export LRGSG_LIB_CBIN="${LRGSG_LIB_CCORE}/bin"
export LRGSG_LIB_CONFIG="${LRGSG_LIB}/config"
export LRGSG_LIB_GT_PATCHES="${LRGSG_LIB}/gt_patches"
export LRGSG_LIB_NX_PATCHES="${LRGSG_LIB}/nx_patches"
export LRGSG_LIB_STOCPROC="${LRGSG_LIB}/stocproc"
export LRGSG_TEST="${LRGSG_ROOT}/test"
export LRGSG_TOOLS="${LRGSG_ROOT}/tools"
export LRGSG_TOOLS_SCRPT="${LRGSG_TOOLS}/bash"
