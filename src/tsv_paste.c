/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#include "cg.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include "tools.h"
#include "tools_bcol.h"
#include "debug.h"
#include "gztools.h"

int copy_line_check(GZFILE *f, int numfields) {
	register int c;
	c = gz_get(f);
	if (c == EOF) {return 0;}
	while (1) {
		if (c == '\n') break;
		if (c == '\t') {
			numfields--;
			if (numfields == 0) {
				fprintf(stderr,"\nfile %s has more columns in a line than header\n",f->filename);
				exit(1);
			}
		}
		putc_unlocked(c,stdout);
		c = gz_get(f);
		if (c == EOF) break;
	}
	while(--numfields) {
		putc_unlocked('\t',stdout);
	}
	return 1;
}

int header(GZFILE *f,int writecomment) {
	int numfields = 1;
	register int c;
	/* read comments */
	while (1) {
		c=gz_get(f);
		if (c != '#') break;
		if (writecomment) putc_unlocked(c,stdout);
		if (c == '\n') continue;
		while ((c=gz_get(f))!=EOF) {
			if (writecomment) putc_unlocked(c,stdout);
			if (c == '\n') break;
		}	
	}
	/* get number of fields */
	while (1) {
		if (c == EOF) break;
		if (c == '\n') break;
		if (c == '\t') numfields++;
		putc_unlocked(c,stdout);
		c=gz_get(f);
	}
	return numfields;
}

int main(int argc, char *argv[]) {
	GZFILE **fa=NULL;
	int *numfieldsa, numfiles;
	register int i;

	if ((argc < 2)) {
		fprintf(stderr,"Format is: tsv_paste file1 file2 ...\n");
		exit(EXIT_FAILURE);
	}
	numfiles = argc -1;
	fa = (GZFILE **)malloc(numfiles*sizeof(FILE *));
	numfieldsa = (int *)malloc(numfiles*sizeof(int));
	for (i = 0 ; i < numfiles ; i++) {
		fa[i] = gz_open(argv[i+1]);
	}
	numfieldsa[0] = header(fa[0],1);
	for (i = 1 ; i < numfiles ; i++) {
		putc_unlocked('\t',stdout);
		numfieldsa[i] = header(fa[i],0);
	}
	putc_unlocked('\n',stdout);
	while(1) {
		if (!copy_line_check(fa[0],numfieldsa[0])) break;
		for (i = 1 ; i < numfiles ; i++) {
			putc_unlocked('\t',stdout);
			copy_line_check(fa[i],numfieldsa[i]);
		}
		putc_unlocked('\n',stdout);
	}
	for (i = 0 ; i < numfiles ; i++) {
		gz_close(fa[i]);
	}
	exit(EXIT_SUCCESS);
}
