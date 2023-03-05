#!/bin/bash
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
dir="$(dirname "$dir")"

echo "Making debug versions of bins in $dir (fast, using tcc)"

cd $dir/src

CC=tcc
COPT="-g -DDEBUG"
CFLAGS="$COPT -D_FILE_OFFSET_BITS=64 -Wall -I$dir/src -I$dir/src/external"
GZTOOLS_C="../src/gztools.c ../src/lz4.c ../src/lz4tools.c ../src/razf.c ../src/bgzf.c ../src/zstdtools.c ../src/external/zstddeclib.c ../src/external/zlib_amalg.c"
GZPOPEN_C=../src/gzpopen.c
BCOL_TOOLS_C=" ../src/tools_bcol.c ../src/xxhash.c"

$CC $CFLAGS -o ../bin/map2besthit ../src/map2besthit.c ../src/tools.c
$CC $CFLAGS -o ../bin/distr2chr ../src/distr2chr.c ../src/hash.c ../src/tools.c
$CC $CFLAGS -o ../bin/map2sv ../src/map2sv.c ../src/tools.c
$CC $CFLAGS -o ../bin/addcols ../src/addcols.c ../src/tools.c
$CC $CFLAGS -o ../bin/reg_subtract ../src/reg_subtract.c ../src/tools.c
$CC $CFLAGS -o ../bin/reg_join ../src/reg_join.c ../src/tools.c
$CC $CFLAGS -o ../bin/covered ../src/covered.c ../src/tools.c
$CC $CFLAGS -o ../bin/getregions ../src/getregions.c ../src/tools.c
$CC $CFLAGS -o ../bin/getregionsbcol ../src/getregionsbcol.c ../src/tools.c
$CC $CFLAGS -o ../bin/getregionsbcol2 ../src/getregionsbcol2.c ../src/tools.c $BCOL_TOOLS_C $GZTOOLS_C
$CC $CFLAGS -o ../bin/reg_annot ../src/reg_annot.c ../src/tools.c
$CC $CFLAGS -o ../bin/var_annot ../src/var_annot.c ../src/tools.c
$CC $CFLAGS -o ../bin/groupby ../src/groupby.c ../src/tools.c
$CC $CFLAGS -o ../bin/multireg ../src/multireg.c ../src/tools.c ../src/tools_var.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/bcol_indexfile ../src/bcol_indexfile.c ../src/tools.c
$CC $CFLAGS -o ../bin/bcol_indexfile_all ../src/bcol_indexfile_all.c ../src/tools.c ../src/hash.c
$CC $CFLAGS -o ../bin/bcol_make ../src/bcol_make.c ../src/tools.c $BCOL_TOOLS_C $GZTOOLS_C
$CC $CFLAGS -o ../bin/bcol_make_multi ../src/bcol_make_multi.c ../src/tools.c $BCOL_TOOLS_C $GZTOOLS_C
$CC $CFLAGS -o ../bin/bcol_annot ../src/bcol_annot.c ../src/tools.c $BCOL_TOOLS_C $GZTOOLS_C
$CC $CFLAGS -o ../bin/test ../src/test.c ../src/tools.c
$CC $CFLAGS -o ../bin/vcf2tsv ../src/vcf2tsv.c ../src/tools.c ../src/hash.c
$CC $CFLAGS -o ../bin/reg_select ../src/reg_select.c ../src/tools.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/depth_histo ../src/depth_histo.c ../src/tools.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/check_sort ../src/check_sort.c ../src/tools.c
$CC $CFLAGS -o ../bin/multi_merge ../src/multi_merge.c ../src/tools.c ../src/tools_var.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/multi_join ../src/multi_join.c ../src/tools.c ../src/tools_var.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/multicompar_addvars ../src/multicompar_addvars.c ../src/tools.c ../src/tools_var.c $BCOL_TOOLS_C $GZTOOLS_C
$CC $CFLAGS -o ../bin/multidb_join ../src/multidb_join.c ../src/tools.c ../src/tools_var.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/tsv_paste ../src/tsv_paste.c ../src/tools.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/tsv_histo ../src/tsv_histo.c ../src/tools.c ../src/hash.c
$CC $CFLAGS -o ../bin/sam_clipamplicons ../src/sam_clipamplicons.c ../src/tools.c
$CC $CFLAGS -o ../bin/tsv2bed ../src/tsv2bed.c ../src/tools.c ../src/tools_var.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/csv2tsv ../src/csv2tsv.c ../src/tools.c
$CC $CFLAGS -o ../bin/tsv2csv ../src/tools.c ../src/tsv2csv.c
$CC $CFLAGS -o ../bin/samcat ../src/samcat.c ../src/tools.c $GZPOPEN_C
$CC $CFLAGS -o ../bin/sam2tsv ../src/sam2tsv.c ../src/tools.c hash.c
$CC $CFLAGS -o ../bin/tsv2sam ../src/tsv2sam.c ../src/tools.c
$CC $CFLAGS -o ../bin/sam_amplicons_count ../src/sam_amplicons_count.c ../src/tools.c
$CC $CFLAGS -o ../bin/varsfile ../src/varsfile.c ../src/tools.c ../src/tools_var.c $GZTOOLS_C
$CC $CFLAGS -o ../bin/splitfastq ../src/splitfastq.c ../src/tools.c
$CC $CFLAGS -o ../bin/countlines ../src/countlines.c ../src/tools.c
$CC $CFLAGS -o ../bin/addrow ../src/addrow.c ../src/tools.c
$CC $CFLAGS -o ../bin/distrreg ../src/distrreg.c ../src/tools.c
$CC $CFLAGS -o ../bin/mergesorted ../src/mergesorted.c ../src/tools.c $GZPOPEN_C
$CC $CFLAGS -o ../bin/test_naturalcompare ../src/test_naturalcompare.c ../src/tools.c
$CC $CFLAGS -o ../bin/noise ../src/tools.c ../src/noise.c
$CC $CFLAGS -o ../bin/estimate_error_rate ../src/tools.c ../src/tools_var.c ../src/estimate_error_rate.c $GZTOOLS_C

echo "Done"
echo "Reminder: $CC included in extra compiles binaries statically linked with musl, so you need the following option to use valgrind:"
echo "valgrind --run-libc-freeres=no bin/exec ..."
