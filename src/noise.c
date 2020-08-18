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
	double numbases=0,nummismatch=0,numins=0,numdel=0;
	int maxcol = 5,count, debug = -1, progress=1000000,progresscur=1000000;
	int mindepth=20,perpos=0,depth=0,nrdiff,i,num;
	int *results;
	if ((argc < 2 || argc > 5)) {
		fprintf(stderr,"Format is: noise mindepth perpos");
		exit(EXIT_FAILURE);
	}
	if (argc >= 2) mindepth=atoi(argv[1]);
	if (argc >= 3) perpos=atoi(argv[2]);
	if (argc >= 4) debug=atoi(argv[3]);
	if (argc >= 5) progress=atoi(argv[4]);
	progresscur = progress;
	results = (int *)malloc((mindepth+2)*sizeof(int));
	for (i = 0 ; i <= mindepth ; i++) {
		results[i] = 0;
	}
	DStringInit(&line);
	DStringGetLine(&line, stdin);
	if (perpos) {
		fprintf(stdout,"depth\tnrdiff\n");
	} else if (debug >= 0) {
		fprintf(stdout,"depth\tnrdiff\tbin\tchromosome\tpos\tref\tdepth\tbases\tqual\n");
	}
	while (1) {
		linepos = line.string;
NODPRINT("%s\n",linepos)
		count = 0;
		if (--progresscur == 0) {
			fprintf(stderr,"%s\n",linepos);
			progresscur = progress;
		}
		while (*linepos && (count <= maxcol)) {
			if (*linepos == '\t') {
				if (*linepos == '\t') {
					scanpos = linepos+1;
				} else {
					scanpos = linepos;
				}
				if (count == 2) {
					depth = atoi(scanpos);
					/* if (depth < mindepth) break; */
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
			while (*linepos != '\t' && *linepos != '\0') {
				char c = *linepos;
				if (c == '+' || c == '-') {
					if (c == '+') {
						nrdiff++;
						numins++;
					} else {
						numdel++;
					}
					linepos++;
					num=atoi(linepos);
					while (*linepos >= '0' && *linepos <= '9') linepos++;
					while (num--) linepos++;
				} else {
					if ((*linepos >= 'A' && *linepos <= 'Z') || (*linepos >= 'a' && *linepos <= 'z')) {
						nrdiff++;
						numbases++; nummismatch++;
					} else if (*linepos == '*') {
						nrdiff++;
					} else if (*linepos == '.' || *linepos == ',') {
						numbases++;
					}
					linepos++;
				}
			}
			if (depth >= mindepth) {
				if (nrdiff > depth) {nrdiff = depth;}
				if (nrdiff == 0) {
					i = 0;
				} else {
					i = (int)(0.9999999999+(mindepth*(double)nrdiff/depth));
				}
				if (i>mindepth) {
					fprintf(stderr,"at line: %s\n",line.string);
					fprintf(stderr,"internal error i > mindepth: i=%d mindepth=%d depth=%d nrdiff=%d\n",i,mindepth,depth,nrdiff);
					exit(1);
				}
				results[i]++;
				if (debug >= 0 && nrdiff >= debug) fprintf(stdout,"%d\t%d\t%.1f\t%s\n",depth,nrdiff,100.0*i/mindepth,line.string);
				if (perpos) fprintf(stdout,"%d\t%d\n",depth,nrdiff);
			}
		}
		read = DStringGetLine(&line, stdin);
		if ((int)read == -1) break;
	}
	if (debug == -1 && perpos == 0) {
		fprintf(stdout,"#filetype\ttsv/noisefile\n");
		fprintf(stdout,"#fileversion\t0.99\n");
		fprintf(stdout,"#description\ta noise (=percentage alternative allele) histogram for all positions with depth >= mindepth \n");
		fprintf(stdout,"#mindepth\t%d\n",mindepth);
		fprintf(stdout,"#general bam stats:\n");
		fprintf(stdout,"#numalignedbases\t%.0f\n",numbases);
		fprintf(stdout,"#nummismatch\t%.0f\n",nummismatch);
		fprintf(stdout,"#numdel\t%.0f\n",numdel);
		fprintf(stdout,"#numins\t%.0f\n",numins);
		fprintf(stdout,"percentage\tcount\n");
		for (i = 0 ; i <= mindepth ; i++) {
			fprintf(stdout,"%0.1f\t%d\n",100.0*i/mindepth,results[i]);
		}
	}
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
