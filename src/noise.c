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
#include "debug.h"

int main(int argc, char *argv[]) {
	DString line;
	char *linepos = NULL, *scanpos = NULL;
	ssize_t read;
	int maxcol = 5,count, debug = -1;
	int mindepth=20,depth,nrdiff,i;
	int *results;
	if ((argc < 2 || argc > 3)) {
		fprintf(stderr,"Format is: noise mindepth");
		exit(EXIT_FAILURE);
	}
	if (argc >= 2) mindepth=atoi(argv[1]);
	if (argc >= 3) debug=atoi(argv[2]);
	results = (int *)malloc((mindepth+1)*sizeof(int));
	for (i = 0 ; i <= mindepth ; i++) {
		results[i] = 0;
	}
	DStringInit(&line);
	DStringGetLine(&line, stdin);
	if (debug >= 0) fprintf(stdout,"depth\tnrdiff\tbracket\tbin\tchromosome\tpos\tref\tdepth\tbases\tqual\n");
	while (1) {
		linepos = line.string;
NODPRINT("%s\n",linepos)
		count = 0;
		while (*linepos && (count <= maxcol)) {
			if (*linepos == '\t') {
				if (*linepos == '\t') {
					scanpos = linepos+1;
				} else {
					scanpos = linepos;
				}
				if (count == 2) {
					depth = atoi(scanpos);
					if (depth < mindepth) break;
				} else if (count == 3) {
					break;
				}
				count++;
			}
			linepos++;
		}
		if (count >= 3) {
			linepos++;
			nrdiff = 0;
			while (*linepos != '\t') {
				char c = *linepos, num;
				if (c == '+' || c == '-') {
					if (c == '+') nrdiff++;
					linepos++;
					num=atoi(linepos);
					while (*linepos >= '0' && *linepos <= '9') linepos++;
					while (num--) linepos++;
				} else {
					if ((*linepos >= 'A' && *linepos <= 'Z') || (*linepos >= 'a' && *linepos <= 'z') || *linepos == '*') nrdiff++;
					linepos++;
				}
			}
			if (nrdiff == 0) {
				i = 0;
			} else {
				i = (int)(0.9999999999+(mindepth*(double)nrdiff/depth));
			}
			results[i]++;
			if (debug >= 0 && nrdiff > debug) fprintf(stdout,"%d\t%d\t%.1f\t%s\n",depth,nrdiff,100.0*i/mindepth,line.string);
		}
		read = DStringGetLine(&line, stdin);
		if ((int)read == -1) break;
	}
	if (debug == -1) {
		fprintf(stdout,"percentage\tcount\n");
		for (i = 0 ; i <= mindepth ; i++) {
			fprintf(stdout,"%0.1f\t%d\n",100.0*i/mindepth,results[i]);
		}
	}
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
