/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tools.h"
#include "gztools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	FILE *f=stdin;
	DString *line;
	size_t read;
	char *cur,prev,c,state;
	if ((argc != 1)) {
		fprintf(stderr,"Format is: csv2tsv\n");
		exit(EXIT_FAILURE);
	}
	/* allocate */
	line = DStringNew();
	prev = '\n';
	state = 'c';
	while (1) {
		read = DStringGetLine(line, f);
		if (read == -1) break;
		cur = line->string;
		while (read--) {
			c = *cur;
			if (c == '"') {
				if (state == 'c') {
					state = 'o';
				} else {
					state = 'c';
				}
				if (prev == '"') {fputc('"',stdout);}
			} else if (c == ',') {
				if (state == 'c') {
					fputc('\t',stdout);
				} else {
					fputc(',',stdout);
				}
			} else if (c == '\t') {
				fputc('\\',stdout);
				fputc('t',stdout);
			} else {
				fputc(c,stdout);
			}
			prev = c;
			cur++;
		}
		if (state == 'c') {
			fputc('\n',stdout);
		} else {
			fputc('\\',stdout);
			fputc('n',stdout);
		}
		prev = '\n';
	}
	if (line) {DStringDestroy(line);}
	exit(EXIT_SUCCESS);
}
