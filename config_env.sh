#!/bin/sh

# Determine the project root directory based on the location of this script
LRGSG_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Export directories
export LRGSG_ROOT
