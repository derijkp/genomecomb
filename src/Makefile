# ARCH should be defined by caller
CFLAGS=-Wall -DDEBUG=1 -g -I~/dev/completegnomics/src
#CFLAGS=-Wall -O2
TOOLS= tools.c tools.h debug.h
HASHTOOLS=hash.h hash.c
BCOL_TOOLS=tools_bcol.c tools_bcol.h razf.c razf.h

all: ../bin$(ARCH)/map2besthit ../bin$(ARCH)/distr2chr ../bin$(ARCH)/map2sv ../bin$(ARCH)/addcols ../bin$(ARCH)/reg_subtract ../bin$(ARCH)/reg_join ../bin$(ARCH)/covered ../bin$(ARCH)/getregions ../bin$(ARCH)/getregionsbcol ../bin$(ARCH)/reg_annot ../bin$(ARCH)/var_annot ../bin$(ARCH)/groupby ../bin$(ARCH)/multireg ../bin$(ARCH)/bcol_indexfile ../bin$(ARCH)/bcol_indexfile_all ../bin$(ARCH)/bcol_make

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

../bin$(ARCH)/reg_annot: reg_annot.c $(TOOLS)
	gcc reg_annot.c tools.c $(CFLAGS) -o ../bin$(ARCH)/reg_annot

../bin$(ARCH)/var_annot: var_annot.c $(TOOLS)
	gcc var_annot.c tools.c $(CFLAGS) -o ../bin$(ARCH)/var_annot

../bin$(ARCH)/groupby: groupby.c $(TOOLS)
	gcc groupby.c tools.c $(CFLAGS) -o ../bin$(ARCH)/groupby

../bin$(ARCH)/multireg: multireg.c $(TOOLS)
	gcc multireg.c tools.c $(CFLAGS) -o ../bin$(ARCH)/multireg

../bin$(ARCH)/bcol_indexfile: bcol_indexfile.c $(TOOLS)
	gcc bcol_indexfile.c tools.c $(CFLAGS) -o ../bin$(ARCH)/bcol_indexfile

../bin$(ARCH)/bcol_indexfile_all: bcol_indexfile_all.c $(TOOLS) $(HASHTOOLS)
	gcc bcol_indexfile_all.c tools.c hash.c $(CFLAGS) -o ../bin$(ARCH)/bcol_indexfile_all

../bin$(ARCH)/bcol_make: bcol_make.c $(TOOLS) $(BCOL_TOOLS) $(HASHTOOLS)
	gcc bcol_make.c tools.c tools_bcol.c hash.c razf.c $(CFLAGS) -o ../bin$(ARCH)/bcol_make -Wl,-Bstatic -lz -Wl,-Bdynamic

../bin$(ARCH)/bcol_annot: bcol_annot.c $(TOOLS) $(BCOL_TOOLS)
	gcc bcol_annot.c tools.c tools_bcol.c razf.c $(CFLAGS) -o ../bin$(ARCH)/bcol_annot -Wl,-Bstatic -lz -Wl,-Bdynamic

../bin$(ARCH)/test: test.c tools.c
	gcc test.c tools.c $(CFLAGS) -o ../bin$(ARCH)/test
