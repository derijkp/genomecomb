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

int main(int argc, char *argv[]) {
	FILE *f = NULL;
	int pos;
	register int c, first = 1;
	if ((argc < 1)) {
		fprintf(stderr,"Format is: samcat samfile1 ...");
		exit(EXIT_FAILURE);
	}
	f = fopen64_or_die(argv[1],"r");
	while ((c=getc_unlocked(f))!=EOF) {
		putc_unlocked(c,stdout);
	}
	fclose(f);
	pos = 2;
	while (pos < argc) {
		f = fopen64_or_die(argv[pos++],"r");
		while ((c=getc_unlocked(f))!=EOF) {
			if (first) {
				if (c != '@') break;
				first = 0;
			}
			if (c == '\n') {first = 1;}
		}
		putc_unlocked(c,stdout);
		while ((c=getc_unlocked(f))!=EOF) {
			putc_unlocked(c,stdout);
		}
		fclose(f);
	}
	exit(EXIT_SUCCESS);
}
