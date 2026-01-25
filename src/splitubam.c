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
	char *output_dir, *output_post, *header;
	uint64_t count;
	unsigned int output_dir_size, output_post_size, numlines,maxparts = INT_MAX, part = 1, threads = 1;
	char c;
	if ((argc < 6)||(argc > 7)) {
		fprintf(stderr,"Format is: splitubam output_dir outputpost numlines header maxparts\n");
		exit(EXIT_FAILURE);
	}
	output_dir = argv[1];
	output_dir_size = strlen(output_dir);
	output_post = argv[2];
	output_post_size = strlen(output_post);
	numlines = atoi(argv[3]);
	header = argv[4];
	threads = atoi(argv[5]);
	if ((argc == 7) && argv[6][0] != '.') {
		maxparts = atoi(argv[6]);
	}
	filename = DStringNew();
	count = 0;
	c = getc_unlocked(stdin);
	while (1) {
		if (c == EOF) break;
		if (count == 0) {
			if (o != NULL) {pclose(o);}
			DStringSetS(filename,"",0);
			DStringPrintf(filename,"samtools view -bhS -@ %d > ",threads);
			fprintf(stdout,"1 %*.*s\n",filename->size,filename->size,filename->string);
			DStringAppendS(filename,output_dir,output_dir_size);
			fprintf(stdout,"2 %*.*s\n",filename->size,filename->size,filename->string);
			DStringPrintf(filename,"/p%d_",part);
			fprintf(stdout,"3 %*.*s\n",filename->size,filename->size,filename->string);
			DStringAppendS(filename,output_post,output_post_size);
			fprintf(stdout,"%*.*s\n",filename->size,filename->size,filename->string);
			o = popen(filename->string,"w");
			fprintf(o,"%s",header);
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
	pclose(o);
	if (filename) {DStringDestroy(filename);}
	exit(EXIT_SUCCESS);
}
