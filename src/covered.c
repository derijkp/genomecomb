/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	FILE *f1;
	DString line;
	DString chromosome1,keepchromosome;
	int64_t size=0,totsize=0;
	int chr1pos,start1pos,end1pos,max1;
	int nchr1=0,start1=-1,end1=-1;
	int error,curchr=0,nextpos=0;
	DStringInit(&line);DStringInit(&chromosome1);DStringInit(&keepchromosome);
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
		error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1);
		if (error) break;
		if (curchr == 0) {
			DStringCopy(&keepchromosome,&chromosome1);
			curchr = 1;
		}
		if (nchr1 > curchr) {
			totsize += size;
			fprintf(stdout,"%s\t%lld\n", keepchromosome.string, size);
			size = 0;
			curchr = nchr1;
			nextpos = 0;
			DStringCopy(&keepchromosome,&chromosome1);
		}
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",chromosome1.string,start1);
			fflush(stderr);
			nextpos += 10000000;
		}
		size += end1 - start1;
	}
	totsize += size;
	fprintf(stdout,"%s\t%lld\n", keepchromosome.string, size);
	fprintf(stdout,"\ntotal\t%lld\n", totsize);
	fclose(f1);
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
