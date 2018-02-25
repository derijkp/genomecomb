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
#include <stdint.h>
#include <inttypes.h>
#include "tools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	FILE *f1;
	DString line;
	DString *chromosome1,keepchromosome;
	int64_t size=0,totsize=0;
	int chr1pos,start1pos,end1pos,max1;
	int start1=-1,end1=-1;
	int error,curchr=0;
#ifdef SHOWPROGRESS
	int nextpos=0;
#endif
	DStringInit(&line);DStringInit(&keepchromosome);
	chromosome1 = DStringNew();
	if ((argc != 4)) {
		fprintf(stderr,"Format is: covered chrpos1 startpos1 endpos1\n");
		exit(EXIT_FAILURE);
	}
	f1 = stdin;
	chr1pos = atoi(argv[1]);
	start1pos = atoi(argv[2]);
	end1pos = atoi(argv[3]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	fprintf(stdout,"chromosome\tbases\n");
	while (1) {
		error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1);
		if (error) break;
		if (curchr == 0) {
			DStringCopy(&keepchromosome,chromosome1);
			curchr = 1;
		}
NODPRINT("%s <-> %s\n",keepchromosome.string,Loc_ChrString(chromosome1))
		if (DStringLocCompare(chromosome1,&keepchromosome)) {
			totsize += size;
			fprintf(stdout,"%s\t%" PRId64 "\n", keepchromosome.string, size);
NODPRINT("is new\n")
			size = 0;
#ifdef SHOWPROGRESS
			nextpos = 0;
#endif
			DStringCopy(&keepchromosome,chromosome1);
		}
#ifdef SHOWPROGRESS
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",Loc_ChrString(chromosome1),start1);
			fflush(stderr);
			nextpos += 10000000;
		}
#endif
		size += end1 - start1;
	}
	totsize += size;
	fprintf(stdout,"%s\t%" PRId64 "\n", keepchromosome.string, size);
	fprintf(stdout,"\ntotal\t%" PRId64 "\n", totsize);
	fclose(f1);
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
