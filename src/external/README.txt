The files in this directory are amalgamations from the zlib and zstd libraries: a condensation of
all source files into a single C and header file, allowing direct incorporation in an exacutable 
during compilation (as e.g. used in the tcc dbg builds)

The zstd amalgamation was made using contrib/single_file_libs/create_single_file_library.sh 
included in the source distribition of zstd v 1.4.8

The zlib amalgamation was made by Forrest Smith (https://www.forrestthewoods.com/blog/improving_open_source_with_amalgamation)
and was downloaded from:
https://github.com/forrestthewoods/lib_fts/blob/master/tests/amalgamate/zlib_amalg.c
https://github.com/forrestthewoods/lib_fts/blob/master/tests/amalgamate/zlib_amalg.h
