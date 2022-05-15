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
	char *linepos = NULL, *scanpos = NULL,*chr=NULL, typediff='a';
	char ref;
	ssize_t read;
	double numbases=0,nummismatch=0,numins=0,numdel=0;
	int maxcol = 5,count, debug = -1, progress=1000000,progresscur=1000000,chrlen=0;
	int mindepth=20,perpos=0,depth=0,nrdiff,nrmismatch,nrmatch,nrins,nrdel,nrA,nrC,nrG,nrT,i,num,begin=0;
	int *results;
	if ((argc < 2 || argc > 6)) {
		fprintf(stderr,"Format is: noise mindepth ?perpos? ?debug? ?progress? ?typediff?");
		exit(EXIT_FAILURE);
	}
	if (argc >= 2) mindepth=atoi(argv[1]);
	if (argc >= 3) perpos=atoi(argv[2]);
	if (argc >= 4) debug=atoi(argv[3]);
	if (argc >= 5) progress=atoi(argv[4]);
	if (argc >= 6) typediff=argv[5][0];
	progresscur = progress;
	results = (int *)malloc((mindepth+2)*sizeof(int));
	for (i = 0 ; i <= mindepth ; i++) {
		results[i] = 0;
	}
	DStringInit(&line);
	DStringGetLine(&line, stdin);
	if (perpos) {
		fprintf(stdout,"chromosome\tbegin\tdepth\tref\tnr_diff\tnr_mismatch\tnr_del\tnr_ins");
		fprintf(stdout,"\tnr_A\tnr_C\tnr_G\tnr_T");
		fprintf(stdout,"\tpct_diff\tpct_mismatch\tpct_del\tpct_ins");
		fprintf(stdout,"\tpct_A\tpct_C\tpct_G\tpct_T\tpct_maxmismatch");
		fprintf(stdout,"\n");
	} else if (debug >= 0) {
		fprintf(stdout,"depth\tnrdiff\tbin\tchromosome\tpos\tref\tdepth\tbases\tqual\n");
	}
	while (1) {
		ref='N';
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
				if (count == 0) {
					chr = line.string;
					chrlen = scanpos - line.string - 1;
					begin = atoi(scanpos)-1;
				} else if (count == 1) {
					ref = *scanpos;
				} else if (count == 2) {
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
			nrdiff = 0; nrmismatch = 0;nrmatch = 0; nrins = 0; nrdel = 0; nrA = 0; nrC = 0 ; nrG = 0; nrT = 0;
			while (*linepos != '\t' && *linepos != '\0') {
				char c = *linepos;
				if (c == '+' || c == '-') {
					if (c == '+') {
						/* insertion after this base */
						nrdiff++; nrins++;
						numins++;
					} else {
						/* next base(s) is a deletion, will be counted in next columns */
					}
					linepos++;
					/* indel symbol (+/-) is followed by:
					  - number of inserted/deleted bases
					  - inserted/deleted bases */
					num=atoi(linepos);
					while (*linepos >= '0' && *linepos <= '9') linepos++;
					while (num--) linepos++;
				} else if (c == '^') {
					/* first position covered by the read , followed by mapping quality encoded as one ascii char*/
					linepos++; linepos++;
				} else {
					if ((*linepos >= 'A' && *linepos <= 'Z') || (*linepos >= 'a' && *linepos <= 'z')) {
						nrdiff++; nrmismatch++;
						numbases++; nummismatch++;
						if (*linepos == 'A' || *linepos == 'a') {
							nrA++;
						} else if (*linepos == 'C' || *linepos == 'c') {
							nrC++;
						} else if (*linepos == 'G' || *linepos == 'g') {
							nrG++;
						} else if (*linepos == 'T' || *linepos == 't' || *linepos == 'U' || *linepos == 'u') {
							nrT++;
						}
					} else if (*linepos == '*' || *linepos == '#') {
						nrdiff++; numdel++; nrdel++;
					} else if (*linepos == '.' || *linepos == ',') {
						numbases++; nrmatch++;
					}
					/* other characters are ignored for counting :
						$: If this is the last position covered by the read
						> and <: Reference skip (due to CIGAR N -> RNA alignment)
					*/
					linepos++;
				}
			}
			if (depth >= mindepth) {
				int diffs, maxmismatch = 0;
				if (perpos == 0) {
					if (typediff == 'a') {
						diffs = nrdiff;
					} else if (typediff == 'm') {
						diffs = nrmismatch;
					} else if (typediff == 'd') {
						diffs = nrdel;
					} else if (typediff == 'i') {
						diffs = nrins;
					} else if (typediff == 'b') {
						diffs = 0;
						if (nrA > diffs) {
							diffs = nrA;
						}
						if (nrC > diffs) {
							diffs = nrC;
						}
						if (nrG > diffs) {
							diffs = nrG;
						}
						if (nrT > diffs) {
							diffs = nrT;
						}
					} else {
						fprintf(stderr,"Wrong parameter for typediff\n");
						exit(1);
					}
					if (diffs > depth) {diffs = depth;}
					if (diffs == 0) {
						i = 0;
					} else {
						i = (int)(0.9999999999+(mindepth*(double)diffs/depth));
					}
					if (i>mindepth) {
						fprintf(stderr,"at line: %s\n",line.string);
						fprintf(stderr,"internal error i > mindepth: i=%d mindepth=%d depth=%d diffs=%d\n",i,mindepth,depth,diffs);
						exit(1);
					}
					results[i]++;
					if (debug >= 0 && nrdiff >= debug) fprintf(stdout,"%d\t%d\t%.2f\t%s\n",depth,nrdiff,100.0*i/mindepth,line.string);
				} else {
					if (nrA > maxmismatch) {
						maxmismatch = nrA;
					}
					if (nrC > maxmismatch) {
						maxmismatch = nrC;
					}
					if (nrG > maxmismatch) {
						maxmismatch = nrG;
					}
					if (nrT > maxmismatch) {
						maxmismatch = nrT;
					}
					if (ref == 'A' || ref == 'a') {
						nrA = nrmatch;
					} else if (ref == 'C' || ref == 'c') {
						nrC = nrmatch;
					} else if (ref == 'G' || ref == 'g') {
						nrC = nrmatch;
					} else if (ref == 'T' || ref == 't' || ref == 'U' || ref == 'u') {
						nrT = nrmatch;
					}
					fprintf(stdout,"%*.*s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
						chrlen,chrlen,chr,begin,depth,ref,
						nrdiff,nrmismatch,nrdel,nrins,
						nrA,nrC,nrG,nrT,
						100.0*nrdiff/depth,100.0*nrmismatch/depth,100.0*nrdel/depth,100.0*nrins/depth,
						100.0*nrA/depth,100.0*nrC/depth,100.0*nrG/depth,100.0*nrT/depth,
						100.0*maxmismatch/depth
					);
				}
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
			fprintf(stdout,"%0.2f\t%d\n",100.0*i/mindepth,results[i]);
		}
	}
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
