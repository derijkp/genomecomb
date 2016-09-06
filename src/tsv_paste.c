/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include "tools.h"
#include "tools_bcol.h"
#include "debug.h"

int copy_line_check(FILE *f, int numfields) {
	register int c;
	c = getc_unlocked(f);
	if (c == EOF) {return 0;}
	while (1) {
		if (c == '\n') break;
		if (c == '\t') {
			numfields--;
			if (numfields == 0) {
				fprintf(stderr,"\nfile has more columns in a line than header\n");
				exit(1);
			}
		}
		putc_unlocked(c,stdout);
		c = getc_unlocked(f);
		if (c == EOF) break;
	}
	while(--numfields) {
		putc_unlocked('\t',stdout);
	}
	return 1;
}

int header(FILE *f,int writecomment) {
	int numfields = 1;
	register int c;
	/* read comments */
	while (1) {
		c=getc_unlocked(f);
		if (c != '#') break;
		putc_unlocked(c,stdout);
		if (c == '\n') continue;
		while ((c=getc_unlocked(f))!=EOF) {
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
		c=getc_unlocked(f);
	}
	return numfields;
}

int main(int argc, char *argv[]) {
	FILE **fa=NULL;
	int *numfieldsa, numfiles;
	register int i;

	if ((argc < 2)) {
		fprintf(stderr,"Format is: tsv_paste file1 file2 ...\n");
		exit(EXIT_FAILURE);
	}
	numfiles = argc -1;
	fa = (FILE **)malloc(numfiles*sizeof(FILE *));
	numfieldsa = (int *)malloc(numfiles*sizeof(int));
	for (i = 0 ; i < numfiles ; i++) {
		fa[i] = fopen64_or_die(argv[i+1],"r");
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
		fclose(fa[i]);
	}
	exit(EXIT_SUCCESS);
}
