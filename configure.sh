
# Source directory of BeatIt. Use full Path
BeatIt_SOURCE_DIR=

# Install directory
PREFIX=

cmake \
  -D CMAKE_ISNTALL_PREFIX=$PREFIX \
  $BeatIt_SOURCE_DIR
