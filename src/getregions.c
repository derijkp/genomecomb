/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	DString line;
	char *chr, *linepos = NULL, *scanpos = NULL;
	ssize_t read;
	int poscol,valuecol,above,shift,maxcol,pos,prev,value,cutoff,accept;
	int status,begin=0,count;
	if ((argc != 7)) {
		fprintf(stderr,"Format is: getregions chromosome poscol valuecol cutoff above shift");
		exit(EXIT_FAILURE);
	}
	DStringInit(&line);
	chr = argv[1];
	poscol = atoi(argv[2]);
	valuecol = atoi(argv[3]);
	cutoff = atoi(argv[4]);
	above = atoi(argv[5]);
	shift = atoi(argv[6]);
	pos = 0 - shift;
	maxcol = poscol ; if (valuecol > maxcol) {maxcol = valuecol;}
	skip_header(stdin,&line);
	DStringGetLine(&line, stdin);
	begin = -1;
	status = 0;
	while (1) {
		linepos = line.string;
NODPRINT("%s\n",linepos)
		count = 0;
		while (*linepos && (count <= maxcol)) {
			if (*linepos == '\t' || (linepos == line.string)) {
				if (*linepos == '\t') {
					scanpos = linepos+1;
				} else {
					scanpos = linepos;
				}
				if (count == poscol) {
					sscanf(scanpos,"%d",&pos);
				} else if (count == valuecol) {
					sscanf(scanpos,"%d",&value);
				}
				count++;
			}
			linepos++;
		}
		if (count > maxcol) {
			accept = ((above && (value > cutoff)) || (!above && (value < cutoff)));
			if (status) {
				if ((prev != pos) || !accept) {
					fprintf(stdout,"%s\t%d\t%d\n", chr, begin+shift, prev+shift);
					status = 0;
					if ((prev != pos) && accept) {
						begin = pos;
						prev = pos+1;
						status = 1;
					}
				} else {
					prev = pos+1;
				}
			} else {
				if (accept) {
					begin = pos;
					prev = pos+1;
					status = 1;
				}
			}
		}
		read = DStringGetLine(&line, stdin);
		if ((int)read == -1) break;
	}
	if (status) {
		fprintf(stdout,"%s\t%d\t%d\n", chr, begin+shift, pos+1+shift);
	}
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
