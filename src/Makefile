# ARCH should be defined by caller
 CFLAGS=-Wall -DDEBUG=1 -g -I~/dev/genomecomb/src
#CFLAGS=-Wall -O2 -I~/dev/genomecomb/src
TOOLS= tools.c tools.h debug.h
GZTOOLS_C= lz4.c lz4tools.c razf.c gztools.c
GZTOOLS= $(GZTOOLS_C) gztools.h razf.h lz4.h lz4tools.h debug.h
TOOLS_VAR= tools_var.c tools_var.h
HASHTOOLS=hash.h hash.c
BCOL_TOOLS=tools_bcol.c tools_bcol.h razf.c razf.h lz4.c lz4.h lz4tools.c lz4tools.h xxhash.c xxhash.h
LIBZ=-Wl,-Bstatic -lz -Wl,-Bdynamic

all: ../bin$(ARCH)/map2besthit ../bin$(ARCH)/distr2chr ../bin$(ARCH)/map2sv ../bin$(ARCH)/addcols \
     ../bin$(ARCH)/reg_subtract ../bin$(ARCH)/reg_join ../bin$(ARCH)/covered ../bin$(ARCH)/getregions \
     ../bin$(ARCH)/getregionsbcol ../bin$(ARCH)/getregionsbcol2 ../bin$(ARCH)/reg_annot ../bin$(ARCH)/var_annot ../bin$(ARCH)/groupby \
     ../bin$(ARCH)/multireg ../bin$(ARCH)/bcol_indexfile ../bin$(ARCH)/bcol_indexfile_all \
     ../bin$(ARCH)/bcol_make ../bin$(ARCH)/bcol_make_multi ../bin$(ARCH)/bcol_annot ../bin$(ARCH)/vcf2tsv ../bin$(ARCH)/reg_select \
     ../bin$(ARCH)/check_sort ../bin$(ARCH)/multi_merge ../bin$(ARCH)/multi_join \
     ../bin$(ARCH)/multicompar_addvars_simple \
     ../bin$(ARCH)/multidb_join \
     ../bin$(ARCH)/tsv_paste ../bin$(ARCH)/tsv_histo \
     ../bin$(ARCH)/sam_clipamplicons

clean:
	rm -f ../bin$(ARCH)/*

../bin$(ARCH)/map2besthit: map2besthit.c $(TOOLS)
	gcc map2besthit.c tools.c $(CFLAGS) -o ../bin$(ARCH)/map2besthit

../bin$(ARCH)/distr2chr: distr2chr.c $(TOOLS) $(HASHTOOLS)
	gcc distr2chr.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/distr2chr

../bin$(ARCH)/map2sv: map2sv.c $(TOOLS)
	gcc map2sv.c tools.c $(CFLAGS) -o ../bin$(ARCH)/map2sv

../bin$(ARCH)/addcols: addcols.c $(TOOLS)
	gcc addcols.c tools.c $(CFLAGS) -o ../bin$(ARCH)/addcols

../bin$(ARCH)/reg_subtract: reg_subtract.c $(TOOLS)
	gcc reg_subtract.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_subtract

../bin$(ARCH)/reg_join: reg_join.c $(TOOLS)
	gcc reg_join.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_join

../bin$(ARCH)/covered: covered.c $(TOOLS)
	gcc covered.c tools.c $(CFLAGS) -o ../bin$(ARCH)/covered

../bin$(ARCH)/getregions: getregions.c $(TOOLS)
	gcc getregions.c tools.c $(CFLAGS) -o ../bin$(ARCH)/getregions

../bin$(ARCH)/getregionsbcol: getregionsbcol.c $(TOOLS)
	gcc getregionsbcol.c tools.c $(CFLAGS) -o ../bin$(ARCH)/getregionsbcol

../bin$(ARCH)/getregionsbcol2: getregionsbcol2.c $(TOOLS)
	gcc getregionsbcol2.c tools.c tools_bcol.c razf.c lz4.c lz4tools.c xxhash.c $(CFLAGS) -o ../bin$(ARCH)/getregionsbcol2 $(LIBZ)

../bin$(ARCH)/reg_annot: reg_annot.c $(TOOLS)
	gcc reg_annot.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_annot

../bin$(ARCH)/var_annot: var_annot.c $(TOOLS)
	gcc var_annot.c tools.c $(CFLAGS) -o ../bin$(ARCH)/var_annot

../bin$(ARCH)/groupby: groupby.c $(TOOLS)
	gcc groupby.c tools.c $(CFLAGS) -o ../bin$(ARCH)/groupby

../bin$(ARCH)/multireg: multireg.c $(TOOLS) $(GZTOOLS) $(TOOLS_VAR)
	gcc multireg.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multireg $(LIBZ)

../bin$(ARCH)/bcol_indexfile: bcol_indexfile.c $(TOOLS)
	gcc bcol_indexfile.c tools.c $(CFLAGS) -o ../bin$(ARCH)/bcol_indexfile

../bin$(ARCH)/bcol_indexfile_all: bcol_indexfile_all.c $(TOOLS) $(HASHTOOLS)
	gcc bcol_indexfile_all.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/bcol_indexfile_all

../bin$(ARCH)/bcol_make: bcol_make.c $(TOOLS) $(BCOL_TOOLS) $(HASHTOOLS)
	gcc bcol_make.c tools.c tools_bcol.c hash.c razf.c lz4.c lz4tools.c xxhash.c $(CFLAGS) -o ../bin$(ARCH)/bcol_make $(LIBZ)

../bin$(ARCH)/bcol_make_multi: bcol_make_multi.c $(TOOLS) $(BCOL_TOOLS) $(HASHTOOLS)
	gcc bcol_make_multi.c tools.c tools_bcol.c hash.c razf.c lz4.c lz4tools.c xxhash.c $(CFLAGS) -o ../bin$(ARCH)/bcol_make_multi $(LIBZ)

../bin$(ARCH)/bcol_annot: bcol_annot.c $(TOOLS) $(BCOL_TOOLS)
	gcc bcol_annot.c tools.c tools_bcol.c razf.c lz4.c lz4tools.c xxhash.c $(CFLAGS) -o ../bin$(ARCH)/bcol_annot $(LIBZ)

../bin$(ARCH)/test: test.c tools.c
	gcc test.c tools.c $(CFLAGS) -o ../bin$(ARCH)/test

../bin$(ARCH)/vcf2tsv: vcf2tsv.c $(TOOLS) $(HASHTOOLS)
	gcc vcf2tsv.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/vcf2tsv

../bin$(ARCH)/reg_select: reg_select.c $(TOOLS)
	gcc reg_select.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_select

../bin$(ARCH)/check_sort: check_sort.c $(TOOLS)
	gcc check_sort.c tools.c $(CFLAGS) -o ../bin$(ARCH)/check_sort

../bin$(ARCH)/multi_merge: multi_merge.c $(TOOLS) $(TOOLS_VAR) $(GZTOOLS)
	gcc multi_merge.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multi_merge $(LIBZ)

../bin$(ARCH)/multi_join: multi_join.c $(TOOLS) $(TOOLS_VAR) $(GZTOOLS)
	gcc multi_join.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multi_join $(LIBZ)

../bin$(ARCH)/multicompar_addvars_simple: multicompar_addvars_simple.c $(TOOLS) $(GZTOOLS) $(TOOLS_VAR)
	gcc multicompar_addvars_simple.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multicompar_addvars_simple $(LIBZ)

../bin$(ARCH)/multidb_join: multidb_join.c $(TOOLS) $(TOOLS_VAR)
	gcc multidb_join.c tools.c tools_var.c $(CFLAGS) -o ../bin$(ARCH)/multidb_join

../bin$(ARCH)/tsv_paste: tsv_paste.c $(TOOLS)
	gcc tsv_paste.c tools.c $(CFLAGS) -o ../bin$(ARCH)/tsv_paste

../bin$(ARCH)/tsv_histo: tsv_histo.c $(TOOLS) $(HASHTOOLS)
	gcc tsv_histo.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/tsv_histo

../bin$(ARCH)/sam_clipamplicons: sam_clipamplicons.c $(TOOLS)
	gcc sam_clipamplicons.c tools.c $(CFLAGS) -o ../bin$(ARCH)/sam_clipamplicons
