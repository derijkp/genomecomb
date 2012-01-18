CFLAGS=-Wall -DDEBUG=1 -g -I~/dev/completegnomics/src
#CFLAGS=-Wall -O2
TOOLS= tools.c tools.h debug.h

all: ../bin/map2besthit ../bin/distr2chr ../bin/map2sv ../bin/addcols ../bin/reg_subtract ../bin/reg_join ../bin/covered ../bin/getregions ../bin/reg_annot ../bin/var_annot ../bin/groupby ../bin/multireg ../bin/bcol_indexfile

clean:
	rm ../bin/*

../bin/map2besthit: map2besthit.c $(TOOLS)
	gcc map2besthit.c tools.c $(CFLAGS) -o ../bin/map2besthit

../bin/distr2chr: distr2chr.c $(TOOLS)
	gcc distr2chr.c tools.c $(CFLAGS) -o ../bin/distr2chr

../bin/map2sv: map2sv.c $(TOOLS)
	gcc map2sv.c tools.c $(CFLAGS) -o ../bin/map2sv

../bin/addcols: addcols.c $(TOOLS)
	gcc addcols.c tools.c $(CFLAGS) -o ../bin/addcols

../bin/reg_subtract: reg_subtract.c $(TOOLS)
	gcc reg_subtract.c tools.c $(CFLAGS) -o ../bin/reg_subtract

../bin/reg_join: reg_join.c $(TOOLS)
	gcc reg_join.c tools.c $(CFLAGS) -o ../bin/reg_join

../bin/covered: covered.c $(TOOLS)
	gcc covered.c tools.c $(CFLAGS) -o ../bin/covered

../bin/getregions: getregions.c $(TOOLS)
	gcc getregions.c tools.c $(CFLAGS) -o ../bin/getregions

../bin/reg_annot: reg_annot.c $(TOOLS)
	gcc reg_annot.c tools.c $(CFLAGS) -o ../bin/reg_annot

../bin/var_annot: var_annot.c $(TOOLS)
	gcc var_annot.c tools.c $(CFLAGS) -o ../bin/var_annot

../bin/groupby: groupby.c $(TOOLS)
	gcc groupby.c tools.c $(CFLAGS) -o ../bin/groupby

../bin/multireg: multireg.c $(TOOLS)
	gcc multireg.c tools.c $(CFLAGS) -o ../bin/multireg

../bin/bcol_indexfile: bcol_indexfile.c $(TOOLS)
	gcc bcol_indexfile.c tools.c $(CFLAGS) -o ../bin/bcol_indexfile

