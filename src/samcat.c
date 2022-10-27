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
#include "tools.h"
#include "gzpopen.h"

int main(int argc, char *argv[]) {
	PFILE *sf = NULL;
	FILE *f = NULL;
	int pos;
	register int c, first = 1;
	if ((argc <= 1)) {
		fprintf(stderr,"Format is: samcat ?-header header? samfile1 ...\n");
		exit(EXIT_FAILURE);
	}
	if (strlen(argv[1]) == 7 && strncmp(argv[1],"-header",7) == 0) {
		f = fopen(argv[2],"r");
		while ((c=getc_unlocked(f))!=EOF) {
			putc_unlocked(c,stdout);
		}
		pos = 3;
	} else {
		sf = gz_popen(argv[1],"sam");
		while ((c=getc_unlocked(sf->f))!=EOF) {
			putc_unlocked(c,stdout);
		}
		gz_pclose(sf);
		pos = 2;
	}
	while (pos < argc) {
		sf = gz_popen(argv[pos++],"sam");
		f = sf->f;
		while ((c=getc_unlocked(f))!=EOF) {
			if (first) {
				if (c != '@') break;
				first = 0;
			}
			if (c == '\n') {first = 1;}
		}
		if (c != EOF) putc_unlocked(c,stdout);
		while ((c=getc_unlocked(f))!=EOF) {
			putc_unlocked(c,stdout);
		}
		gz_pclose(sf);
	}
	exit(EXIT_SUCCESS);
}
