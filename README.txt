MARCOS modular source
=====================

Directory layout:
  main.cpp       Main program / run orchestration
  include/       Shared headers
  src/           C++ implementation files

Time ordering:
  For every run, CoincidenceProcessor first calls sort_module_file() for each
  module input. Each module is sorted independently using bounded-memory
  external sorting (64 MiB chunks by default). Sorted chunks are merged into
  one sorted temporary file per module. CoincidenceProcessor then performs a
  k-way merge of those sorted module streams using one look-ahead event per
  module before constructing coincidence blocks.

  Temporary sort files are written to the system temporary directory and
  removed after processing. Sorting needs temporary disk space approximately
  equal to the input data size, in addition to space for intermediate chunks.
  If any module fails sorting, processing of that run is aborted to avoid
  silently producing incomplete output.

Build:
  make clean
  make -j$(nproc)

Debug build:
  make debug

The build requires ROOT (root-config) and OpenMP. This implementation was
syntax-checked and the external sorter was tested with a small synthetic
unsorted PIXIE-style file. A full ROOT-linked build and comparison against
production data should still be performed locally.
