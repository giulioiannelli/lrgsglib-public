#!/bin/sh

# Array of path variables
paths=(
    "LRGSG_ROOT"
    "LRGSG_BUILD=\"\${LRGSG_ROOT}/build\""
    "LRGSG_LOG=\"\${LRGSG_ROOT}/log\""
    "LRGSG_DATA=\"\${LRGSG_ROOT}/data\""
    "LRGSG_IPYNB=\"\${LRGSG_ROOT}/ipynb\""
    "LRGSG_SRC=\"\${LRGSG_ROOT}/src\""
    "LRGSG_LIB=\"\${LRGSG_SRC}/lrgsglib\""
    "LRGSG_LIB_CCORE=\"\${LRGSG_LIB}/Ccore\""
    "LRGSG_LIB_CBIN=\"\${LRGSG_LIB_CCORE}/bin\""
    "LRGSG_LIB_CONFIG=\"\${LRGSG_LIB}/config\""
    "LRGSG_LIB_GT_PATCHES=\"\${LRGSG_LIB}/gt_patches\""
    "LRGSG_LIB_NX_PATCHES=\"\${LRGSG_LIB}/nx_patches\""
    "LRGSG_LIB_STOCPROC=\"\${LRGSG_LIB}/stocproc\""
    "LRGSG_TEST=\"\${LRGSG_ROOT}/test\""
    "LRGSG_TOOLS=\"\${LRGSG_ROOT}/tools\""
    "LRGSG_TOOLS_SCRPT=\"\${LRGSG_TOOLS}/bash\""
)
