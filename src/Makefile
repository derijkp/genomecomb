CFLAGS=-Wall -g
#CFLAGS=-Wall -O2

all: ../bin/map2besthit ../bin/distr2chr ../bin/map2sv ../bin/addcols ../bin/reg_subtract ../bin/reg_join ../bin/covered ../bin/getregions

clean:
	rm ../bin/*

../bin/map2besthit: map2besthit.c
	gcc map2besthit.c $(CFLAGS) -o ../bin/map2besthit

../bin/distr2chr: distr2chr.c
	gcc distr2chr.c $(CFLAGS) -o ../bin/distr2chr

../bin/map2sv: map2sv.c
	gcc map2sv.c $(CFLAGS) -o ../bin/map2sv

../bin/addcols: addcols.c
	gcc addcols.c $(CFLAGS) -o ../bin/addcols

../bin/reg_subtract: reg_subtract.c
	gcc reg_subtract.c $(CFLAGS) -o ../bin/reg_subtract

../bin/reg_join: reg_join.c
	gcc reg_join.c $(CFLAGS) -o ../bin/reg_join

../bin/covered: covered.c
	gcc covered.c $(CFLAGS) -o ../bin/covered

../bin/getregions: getregions.c
	gcc getregions.c $(CFLAGS) -o ../bin/getregions
