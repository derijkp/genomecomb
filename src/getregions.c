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
	char curchr[1000], *linepos = NULL, *scanpos = NULL;
	ssize_t read;
	int chrcol,poscol,valuecol,above,shift,maxcol,pos,prev=0,value,cutoff,accept,header;
	int status,begin=0,count;
	if ((argc != 9)) {
		fprintf(stderr,"Format is: getregions chromosome chrpos poscol valuecol cutoff above shift header");
		exit(EXIT_FAILURE);
	}
	DStringInit(&line);
	strncpy(curchr,argv[1],999);
	chrcol = atoi(argv[2]);
	poscol = atoi(argv[3]);
	valuecol = atoi(argv[4]);
	cutoff = atoi(argv[5]);
	above = atoi(argv[6]);
	shift = atoi(argv[7]);
	header = atoi(argv[8]);
	pos = 0 - shift;
	maxcol = poscol;
	if (valuecol > maxcol) {maxcol = valuecol;}
	if (chrcol > maxcol) {maxcol = chrcol;}
	if (header) skip_header(stdin,&line,NULL,NULL);
	DStringGetLine(&line, stdin);
	begin = -1;
	/* if status is 1, we are in a region */
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
				} else if (count == chrcol) {
					char *testpos = linepos, *curpos = curchr; int count = 999;
					while (*testpos && *testpos != '\t') {
						if (*testpos != *curpos) {
							if (status) {
								fprintf(stdout,"%s\t%d\t%d\n", curchr, begin+shift, pos+1+shift);
								status = 0;
							}
							while (*testpos && *testpos != '\t') {
								*curpos++ = *testpos++; count--;
								if (!count) break;
							}
							*curpos = '\0';
							break;
						}
						if (*testpos == '\0') break;
						testpos++; curpos++; count--;
						if (!count) break;
					}
				}
				count++;
			}
			linepos++;
		}
		if (count > maxcol) {
			accept = ((above && (value > cutoff)) || (!above && (value < cutoff)));
			if (status) {
				if ((prev != pos) || !accept) {
					fprintf(stdout,"%s\t%d\t%d\n", curchr, begin+shift, prev+shift);
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
		fprintf(stdout,"%s\t%d\t%d\n", curchr, begin+shift, pos+1+shift);
	}
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
