# ARCH should be defined by caller
ARCH ?=
arch ?= x86_64
COPT ?= -O2
GCC ?= gcc
ZSTDDIR ?= $(HOME)/build/bin-$(arch)/zstd-1.4.4
EXEEXT ?=
CFLAGS := -D_FILE_OFFSET_BITS=64 -Wall $(COPT) -I~/dev/genomecomb/src ${CFLAGS} -I$(ZSTDDIR)/lib -L$(ZSTDDIR)/lib
TOOLS = tools.c tools.h debug.h
GZTOOLS_C = lz4.c lz4tools.c razf.c gztools.c zstdtools.c
GZTOOLS = $(GZTOOLS_C) gztools.h razf.h lz4.h lz4tools.h debug.h zstdtools.h
GZPOPEN_C = gzpopen.c
GZPOPEN = $(GZPOPEN_C) gzpopen.h
TOOLS_VAR = tools_var.c tools_var.h
HASHTOOLS = hash.h hash.c
BCOL_TOOLS_C = tools_bcol.c xxhash.c
BCOL_TOOLS = $(BCOL_TOOLS_C) tools_bcol.h xxhash.h
LIBZ=-Wl,-Bstatic -lz -lzstd -Wl,-Bdynamic

all: ../bin$(ARCH)/map2besthit$(EXEEXT) ../bin$(ARCH)/distr2chr$(EXEEXT) ../bin$(ARCH)/map2sv$(EXEEXT) ../bin$(ARCH)/addcols$(EXEEXT) \
     ../bin$(ARCH)/reg_subtract$(EXEEXT) ../bin$(ARCH)/reg_join$(EXEEXT) ../bin$(ARCH)/covered$(EXEEXT) ../bin$(ARCH)/getregions$(EXEEXT) \
     ../bin$(ARCH)/getregionsbcol$(EXEEXT) ../bin$(ARCH)/getregionsbcol2$(EXEEXT) ../bin$(ARCH)/reg_annot$(EXEEXT) ../bin$(ARCH)/var_annot$(EXEEXT) ../bin$(ARCH)/groupby$(EXEEXT) \
     ../bin$(ARCH)/multireg$(EXEEXT) ../bin$(ARCH)/bcol_indexfile$(EXEEXT) ../bin$(ARCH)/bcol_indexfile_all$(EXEEXT) \
     ../bin$(ARCH)/bcol_make$(EXEEXT) ../bin$(ARCH)/bcol_make_multi$(EXEEXT) ../bin$(ARCH)/bcol_annot$(EXEEXT) ../bin$(ARCH)/vcf2tsv$(EXEEXT) \
     ../bin$(ARCH)/reg_select$(EXEEXT) ../bin$(ARCH)/depth_histo$(EXEEXT) \
     ../bin$(ARCH)/check_sort$(EXEEXT) ../bin$(ARCH)/multi_merge$(EXEEXT) ../bin$(ARCH)/multi_join$(EXEEXT) \
     ../bin$(ARCH)/multicompar_addvars$(EXEEXT) \
     ../bin$(ARCH)/multidb_join$(EXEEXT) \
     ../bin$(ARCH)/tsv_paste$(EXEEXT) ../bin$(ARCH)/tsv_histo$(EXEEXT) \
     ../bin$(ARCH)/sam_clipamplicons$(EXEEXT) \
     ../bin$(ARCH)/tsv2bed$(EXEEXT) \
     ../bin$(ARCH)/samcat$(EXEEXT) \
     ../bin$(ARCH)/sam_amplicons_count$(EXEEXT) \
     ../bin$(ARCH)/varsfile$(EXEEXT) \
     ../bin$(ARCH)/sam2tsv$(EXEEXT) ../bin$(ARCH)/tsv2sam$(EXEEXT) \
     ../bin$(ARCH)/splitfastq$(EXEEXT) \
     ../bin$(ARCH)/countlines$(EXEEXT) \
     ../bin$(ARCH)/addrow$(EXEEXT) \
     ../bin$(ARCH)/csv2tsv$(EXEEXT) \
     ../bin$(ARCH)/distrreg$(EXEEXT) \
     ../bin$(ARCH)/mergesorted$(EXEEXT) \
     ../bin$(ARCH)/noise$(EXEEXT)

clean:
	rm -f ../bin$(ARCH)/*

../bin$(ARCH)/map2besthit$(EXEEXT): map2besthit.c $(TOOLS)
	$(GCC) map2besthit.c tools.c $(CFLAGS) -o ../bin$(ARCH)/map2besthit$(EXEEXT)

../bin$(ARCH)/distr2chr$(EXEEXT): distr2chr.c $(TOOLS) $(HASHTOOLS)
	$(GCC) distr2chr.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/distr2chr$(EXEEXT)

../bin$(ARCH)/map2sv$(EXEEXT): map2sv.c $(TOOLS)
	$(GCC) map2sv.c tools.c $(CFLAGS) -o ../bin$(ARCH)/map2sv$(EXEEXT)

../bin$(ARCH)/addcols$(EXEEXT): addcols.c $(TOOLS)
	$(GCC) addcols.c tools.c $(CFLAGS) -o ../bin$(ARCH)/addcols$(EXEEXT)

../bin$(ARCH)/reg_subtract$(EXEEXT): reg_subtract.c $(TOOLS)
	$(GCC) reg_subtract.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_subtract$(EXEEXT)

../bin$(ARCH)/reg_join$(EXEEXT): reg_join.c $(TOOLS)
	$(GCC) reg_join.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_join$(EXEEXT)

../bin$(ARCH)/covered$(EXEEXT): covered.c $(TOOLS)
	$(GCC) covered.c tools.c $(CFLAGS) -o ../bin$(ARCH)/covered$(EXEEXT)

../bin$(ARCH)/getregions$(EXEEXT): getregions.c $(TOOLS)
	$(GCC) getregions.c tools.c $(CFLAGS) -o ../bin$(ARCH)/getregions$(EXEEXT)

# getregionsbcol one does not use/need bcoltools/compression
../bin$(ARCH)/getregionsbcol$(EXEEXT): getregionsbcol.c $(TOOLS)
	$(GCC) getregionsbcol.c tools.c $(CFLAGS) -o ../bin$(ARCH)/getregionsbcol$(EXEEXT)

../bin$(ARCH)/getregionsbcol2$(EXEEXT): getregionsbcol2.c $(TOOLS) $(BCOL_TOOLS) $(GZTOOLS)
	$(GCC) getregionsbcol2.c tools.c $(BCOL_TOOLS_C) $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/getregionsbcol2$(EXEEXT) $(LIBZ)

../bin$(ARCH)/reg_annot$(EXEEXT): reg_annot.c $(TOOLS)
	$(GCC) reg_annot.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_annot$(EXEEXT)

../bin$(ARCH)/var_annot$(EXEEXT): var_annot.c $(TOOLS)
	$(GCC) var_annot.c tools.c $(CFLAGS) -o ../bin$(ARCH)/var_annot$(EXEEXT)

../bin$(ARCH)/groupby$(EXEEXT): groupby.c $(TOOLS)
	$(GCC) groupby.c tools.c $(CFLAGS) -o ../bin$(ARCH)/groupby$(EXEEXT)

../bin$(ARCH)/multireg$(EXEEXT): multireg.c $(TOOLS) $(GZTOOLS) $(TOOLS_VAR)
	$(GCC) multireg.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multireg$(EXEEXT) $(LIBZ)

../bin$(ARCH)/bcol_indexfile$(EXEEXT): bcol_indexfile.c $(TOOLS)
	$(GCC) bcol_indexfile.c tools.c $(CFLAGS) -o ../bin$(ARCH)/bcol_indexfile$(EXEEXT)

../bin$(ARCH)/bcol_indexfile_all$(EXEEXT): bcol_indexfile_all.c $(TOOLS) $(HASHTOOLS)
	$(GCC) bcol_indexfile_all.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/bcol_indexfile_all$(EXEEXT)

../bin$(ARCH)/bcol_make$(EXEEXT): bcol_make.c $(TOOLS) $(BCOL_TOOLS) $(GZTOOLS) $(HASHTOOLS)
	$(GCC) bcol_make.c tools.c $(BCOL_TOOLS_C) $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/bcol_make$(EXEEXT) $(LIBZ)

../bin$(ARCH)/bcol_make_multi$(EXEEXT): bcol_make_multi.c $(TOOLS) $(BCOL_TOOLS) $(GZTOOLS) $(HASHTOOLS)
	$(GCC) bcol_make_multi.c tools.c $(BCOL_TOOLS_C) $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/bcol_make_multi$(EXEEXT) $(LIBZ)

../bin$(ARCH)/bcol_annot$(EXEEXT): bcol_annot.c $(TOOLS) $(BCOL_TOOLS) $(GZTOOLS)
	$(GCC) bcol_annot.c tools.c $(BCOL_TOOLS_C) $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/bcol_annot$(EXEEXT) $(LIBZ)

../bin$(ARCH)/test$(EXEEXT): test.c tools.c
	$(GCC) test.c tools.c $(CFLAGS) -o ../bin$(ARCH)/test$(EXEEXT)

../bin$(ARCH)/vcf2tsv$(EXEEXT): vcf2tsv.c $(TOOLS) $(HASHTOOLS)
	$(GCC) vcf2tsv.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/vcf2tsv$(EXEEXT)

../bin$(ARCH)/reg_select$(EXEEXT): reg_select.c $(TOOLS) $(GZTOOLS)
	$(GCC) reg_select.c tools.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/reg_select$(EXEEXT) $(LIBZ)

../bin$(ARCH)/depth_histo$(EXEEXT): depth_histo.c $(TOOLS) $(GZTOOLS)
	$(GCC) depth_histo.c tools.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/depth_histo$(EXEEXT) $(LIBZ)

../bin$(ARCH)/check_sort$(EXEEXT): check_sort.c $(TOOLS)
	$(GCC) check_sort.c tools.c $(CFLAGS) -o ../bin$(ARCH)/check_sort$(EXEEXT)

../bin$(ARCH)/multi_merge$(EXEEXT): multi_merge.c $(TOOLS) $(TOOLS_VAR) $(GZTOOLS)
	$(GCC) multi_merge.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multi_merge$(EXEEXT) $(LIBZ)

../bin$(ARCH)/multi_join$(EXEEXT): multi_join.c $(TOOLS) $(TOOLS_VAR) $(GZTOOLS)
	$(GCC) multi_join.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multi_join$(EXEEXT) $(LIBZ)

../bin$(ARCH)/multicompar_addvars$(EXEEXT): multicompar_addvars.c $(TOOLS) $(BCOL_TOOLS) $(GZTOOLS) $(TOOLS_VAR)
	$(GCC) multicompar_addvars.c tools.c tools_var.c $(BCOL_TOOLS_C) $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multicompar_addvars$(EXEEXT) $(LIBZ)

../bin$(ARCH)/multidb_join$(EXEEXT): multidb_join.c $(TOOLS) $(TOOLS_VAR) $(GZTOOLS)
	$(GCC) multidb_join.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/multidb_join$(EXEEXT) $(LIBZ)

../bin$(ARCH)/tsv_paste$(EXEEXT): tsv_paste.c $(TOOLS) $(GZTOOLS)
	$(GCC) tsv_paste.c tools.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/tsv_paste$(EXEEXT) $(LIBZ)

../bin$(ARCH)/tsv_histo$(EXEEXT): tsv_histo.c $(TOOLS) $(HASHTOOLS)
	$(GCC) tsv_histo.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/tsv_histo$(EXEEXT)

../bin$(ARCH)/sam_clipamplicons$(EXEEXT): sam_clipamplicons.c $(TOOLS)
	$(GCC) sam_clipamplicons.c tools.c $(CFLAGS) -o ../bin$(ARCH)/sam_clipamplicons$(EXEEXT)

../bin$(ARCH)/tsv2bed$(EXEEXT): tsv2bed.c $(TOOLS) $(TOOLS_VAR) $(GZTOOLS)
	$(GCC) tsv2bed.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/tsv2bed$(EXEEXT) $(LIBZ)

../bin$(ARCH)/csv2tsv$(EXEEXT): csv2tsv.c $(TOOLS)
	$(GCC) csv2tsv.c tools.c $(CFLAGS) -o ../bin$(ARCH)/csv2tsv$(EXEEXT)

../bin$(ARCH)/samcat$(EXEEXT): samcat.c $(TOOLS) $(GZPOPEN)
	$(GCC) samcat.c tools.c $(GZPOPEN_C) $(CFLAGS) -o ../bin$(ARCH)/samcat$(EXEEXT)

../bin$(ARCH)/sam2tsv$(EXEEXT): sam2tsv.c $(TOOLS) $(HASHTOOLS)
	$(GCC) sam2tsv.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/sam2tsv$(EXEEXT)

../bin$(ARCH)/tsv2sam$(EXEEXT): tsv2sam.c $(TOOLS)
	$(GCC) tsv2sam.c tools.c $(CFLAGS) -o ../bin$(ARCH)/tsv2sam$(EXEEXT)

../bin$(ARCH)/sam_amplicons_count$(EXEEXT): sam_amplicons_count.c $(TOOLS)
	$(GCC) sam_amplicons_count.c tools.c $(CFLAGS) -o ../bin$(ARCH)/sam_amplicons_count$(EXEEXT)

../bin$(ARCH)/varsfile$(EXEEXT): varsfile.c $(TOOLS) $(TOOLS_VAR) $(GZTOOLS)
	$(GCC) varsfile.c tools.c tools_var.c $(GZTOOLS_C) $(CFLAGS) -o ../bin$(ARCH)/varsfile$(EXEEXT) $(LIBZ)

../bin$(ARCH)/splitfastq$(EXEEXT): splitfastq.c $(TOOLS)
	$(GCC) splitfastq.c tools.c $(CFLAGS) -o ../bin$(ARCH)/splitfastq$(EXEEXT)

../bin$(ARCH)/countlines$(EXEEXT): countlines.c $(TOOLS)
	$(GCC) countlines.c tools.c $(CFLAGS) -o ../bin$(ARCH)/countlines$(EXEEXT)

../bin$(ARCH)/addrow$(EXEEXT): addrow.c $(TOOLS)
	$(GCC) addrow.c tools.c $(CFLAGS) -o ../bin$(ARCH)/addrow$(EXEEXT)

../bin$(ARCH)/distrreg$(EXEEXT): distrreg.c $(TOOLS)
	$(GCC) distrreg.c tools.c $(CFLAGS) -o ../bin$(ARCH)/distrreg$(EXEEXT)

../bin$(ARCH)/mergesorted$(EXEEXT): mergesorted.c $(TOOLS) $(GZPOPEN)
	$(GCC) mergesorted.c tools.c $(GZPOPEN_C) $(CFLAGS) -o ../bin$(ARCH)/mergesorted$(EXEEXT)

../bin$(ARCH)/test_naturalcompare$(EXEEXT): test_naturalcompare.c tools.c
	echo "test test_naturalcompare"
	$(GCC) test_naturalcompare.c tools.c $(CFLAGS) -o ../bin$(ARCH)/test_naturalcompare$(EXEEXT)

../bin$(ARCH)/noise$(EXEEXT): noise.c $(TOOLS)
	$(GCC) noise.c tools.c $(CFLAGS) -o ../bin$(ARCH)/noise$(EXEEXT)

