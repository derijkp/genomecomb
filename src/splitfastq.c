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
#include <limits.h>
#include "tools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	FILE *o = NULL;
	DString *filename = NULL;
	char *output_dir, *output_post;
	uint64_t count;
	unsigned int output_dir_size, output_post_size, numlines,maxparts = INT_MAX, part = 1;
	char c;
	if ((argc < 4)||(argc > 5)) {
		fprintf(stderr,"Format is: splitfastq output_dir outputpost numlines maxparts\n");
		exit(EXIT_FAILURE);
	}
	output_dir = argv[1];
	output_dir_size = strlen(output_dir);
	output_post = argv[2];
	output_post_size = strlen(output_post);
	numlines = 4*atoi(argv[3]);
	if ((argc == 5) && argv[4][0] != '.') {
		maxparts = atoi(argv[4]);
	}
	filename = DStringNew();
	count = 0;
	c = getc_unlocked(stdin);
	while (1) {
		if (c == EOF) break;
		if (count == 0) {
			if (o != NULL) {fclose(o);}
			DStringSetS(filename,output_dir,output_dir_size);
			DStringPrintf(filename,"/p%d_",part);
			DStringAppendS(filename,output_post,output_post_size);
			fprintf(stdout,"%*.*s\n",filename->size,filename->size,filename->string);
			o = fopen(filename->string,"w");
			if (part < maxparts) {
				count = numlines;
			} else {
				count = UINT64_MAX;
			}
			part++;
		}
		if (c == '\n') count--;
		putc_unlocked(c,o);
		c = getc_unlocked(stdin);
	}
	fclose(o);
	if (filename) {DStringDestroy(filename);}
	exit(EXIT_SUCCESS);
}
