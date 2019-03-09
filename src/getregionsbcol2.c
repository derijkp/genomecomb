/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "tools.h"
#include "tools_bcol.h"
#include "debug.h"


int main(int argc, char *argv[]) {
	BCol *bcol;
	BCol_table *table;
	DString *chr = NULL;
	long double value=0,cutoff=0;
	uint64_t togo;
	int tablepos = 0,tablesize;
	int above=0,shift,pos=0,accept;
	int status,begin=0,error;
	if ((argc != 5)) {
		fprintf(stderr,"Format is: getregions bcolfile cutoff above shift");
		exit(EXIT_FAILURE);
	}
	bcol = bcol_open(argv[1]);
	table = bcol->table;
	tablesize = bcol->tablesize;
	cutoff = atoll(argv[2]);
	above = atoi(argv[3]);
	shift = atoi(argv[4]);
	begin = -1;
	status = 0;
	tablepos = 0; togo = 0;
	while (1) {
		if (!togo) {
			if (status) {
				fprintf(stdout,"%*.*s\t%d\t%d\n", chr->size, chr->size, chr->string, begin, pos);
				status = 0;
			}
			if (tablepos >= tablesize) break;
			chr = table[tablepos].chr;
			togo = table[tablepos].end - table[tablepos].begin;
			pos = table[tablepos].begin + shift;
			tablepos++;
		}
		error = bcol_readdouble(bcol,&value);
		NODPRINT("%d -> %Lf",pos,value)
		if (error == 0) break;
		accept = ((above && (value > cutoff)) || (!above && (value < cutoff)));
		if (status) {
			if (!accept) {
				fprintf(stdout,"%*.*s\t%d\t%d\n", chr->size, chr->size, chr->string, begin, pos);
				status = 0;
			}
		} else {
			if (accept) {
				begin = pos;
				status = 1;
			}
		}
		pos++;
		togo--;
	}
	if (status) {
		fprintf(stdout,"%*.*s\t%d\t%d\n", chr->size, chr->size, chr->string, begin, pos);
	}
	exit(EXIT_SUCCESS);
}
