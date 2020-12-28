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
#include "debug.h"

void write_field(char *prev,char *cur, int quotes) {
	if (quotes) fputc('\"',stdout);
	while (prev < cur) {
		fputc(*prev++,stdout);
	}
	if (quotes) fputc('\"',stdout);
	quotes = 0;
}

int main(int argc, char *argv[]) {
	FILE *f=stdin;
	DString *line;
	int nrread,quotes;
	char *cur,*prev,c;
	if ((argc != 1)) {
		fprintf(stderr,"Format is: tsv2csv\n");
		exit(EXIT_FAILURE);
	}
	/* allocate */
	line = DStringNew();
	while (1) {
		nrread = DStringGetLine(line, f);
		if (nrread == -1) break;
		cur = line->string;
		prev = cur;
		quotes = 0;
		while (nrread--) {
			c = *cur;
			if (c == ',') {
				quotes = 1;
				cur++;
			} else if (c == '\t') {
				write_field(prev,cur, quotes);
				cur++;
				prev=cur;
				fputc(',',stdout);
				quotes = 0;
			} else if (c == '\n') {
				write_field(prev,cur, quotes);
				cur++;
				prev=cur;
				fputc('\n',stdout);
				quotes = 0;
			} else {
				cur++;
			}
		}
		write_field(prev,cur, quotes);
		fputc('\n',stdout);
	}
	if (line) {DStringDestroy(line);}
	exit(EXIT_SUCCESS);
}
