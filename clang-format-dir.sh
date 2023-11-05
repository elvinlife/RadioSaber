#!/bin/bash

# In-place clang-format all .cpp and .h files under a specified directory.
# By default formats all .cpp and .h files under ./
#
# Usage ./clang-format-dir.sh [directory] [--verbose]
#
# Example usage:
#   ./clang-format-dir.sh                   Format all under ./
#   ./clang-format-dir.sh ./src/            Format all under ./src/
#   ./clang-format-dir.sh ./src/ --verbose  Format all under ./src/ and print verbose

directory=${1:-./}
verbose=""

check_directory() {
  local dir=$1
  if [ ! -d "$dir" ]; then
    echo "Directory '$dir' does not exist."
    exit 1
  fi
}

usage() {
  echo "Invalid argument. Usage $0 [directory] [--verbose]"
  exit 1
}

if [ "$#" -eq 1 ]; then
  if [ "$1" == "--verbose" ]; then
    verbose="--verbose"
  else
    directory="$1"
    check_directory "$directory"
  fi
elif [ "$#" -eq 2 ]; then
  if [ "$2" == "--verbose" ]; then
    verbose="--verbose"
    directory="$1"
    check_directory "$directory"
  else
    usage
  fi
else
  usage
fi

find "$directory" -iname "*.h" -o -iname "*.cpp" | xargs clang-format $verbose -i --style=file
