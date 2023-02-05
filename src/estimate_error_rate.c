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
#include "tools_var.h"
#include "gztools.h"
#include "debug.h"

struct sum_info {
	int quality;
	double nummatch;
	double nummismatch;
	double numdel;
	double numins;
};

int main(int argc, char *argv[]) {
	VarFile *varfile;
	DString line;
	char *linepos = NULL, *bases = NULL, *qs = NULL, *scanpos = NULL,*chr=NULL;
	struct sum_info sum_info[60];
	ssize_t read;
	double numbases=0,nummatch=0,nummismatch=0,numins=0,numdel=0;
	int comp;
	int q;
	int maxcol = 5,count, debug = 0, progress=1000000,progresscur=1000000,chrlen=0;
	int nrdiff,nrmismatch,nrmatch,nrins,nrdel,nrA,nrC,nrG,nrT,i,num,begin=0;
	if ((argc < 2 || argc > 4)) {
		fprintf(stderr,"Format is: estimate_error_rate regfile ?debug? ?progress?\n");
		exit(EXIT_FAILURE);
	}
	varfile = OpenVarfile(argv[1],1);
	varfile_next(varfile,1);
	if (argc >= 3) debug=atoi(argv[2]);
	if (argc >= 4) progress=atoi(argv[3]);
	progresscur = progress;
	for ( i = 0 ; i < 60 ; i++) {
		sum_info[i].nummatch = 0;
		sum_info[i].nummismatch = 0;
		sum_info[i].numdel = 0;
		sum_info[i].numins = 0;
	}
	DStringInit(&line);
	DStringGetLine(&line, stdin);
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
				if (count == 0) {
					chr = line.string;
					chrlen = scanpos - line.string - 1;
					begin = atoi(scanpos)-1;
				} else if (count == 1) {
				} else if (count == 2) {
				} else if (count == 3) {
					bases=linepos;
				} else if (count == 4) {
					qs=linepos;
					break;
				}
				count++;
			}
			linepos++;
		}

		while (1) {
			if (varfile->error) break;
			comp = loccompare(chr,varfile->var->chr->string,chrlen,varfile->var->chr->size);
			if (comp == 0) {
				if (begin < varfile->var->start) {
					comp = -1;
				} else if (begin >= varfile->var->end) {
					comp = 1;
				}
			}
			if (comp <= 0) break;
			varfile_next(varfile,1);	
		}
		if (varfile->error) {
			break;
		}
		if (comp > 0) continue;
		if (comp == 0 && count >= 4) {
			if (debug) {
				fprintf(stderr,"%s\t%d\t%s\n",chr,begin,bases);
			}
			bases++;
			qs++;
			nrdiff = 0; nrmismatch = 0;nrmatch = 0; nrins = 0; nrdel = 0; nrA = 0; nrC = 0 ; nrG = 0; nrT = 0;
			while (*bases != '\t' && *bases != '\0') {
				char c = *bases;
				int q = (int)*qs - 33;
				if (c == '+' || c == '-') {
					if (c == '+') {
						/* insertion after this base */
						sum_info[q].numins++;
						nrdiff++; nrins++;
						numins++;
					} else {
						sum_info[q].numdel++;
						numdel++;
					}
					bases++;
					/* indel symbol (+/-) is followed by:
					  - number of inserted/deleted bases
					  - inserted/deleted bases */
					num=atoi(bases);
					while (*bases >= '0' && *bases <= '9') bases++;
					while (num--) bases++;
				} else if (c == '^') {
					/* first position covered by the read , followed by mapping quality encoded as one ascii char*/
					bases++; bases++;
				} else {
					if ((*bases >= 'A' && *bases <= 'Z') || (*bases >= 'a' && *bases <= 'z')) {
						sum_info[q].nummismatch++;
						nrdiff++; nrmismatch++;
						numbases++; nummismatch++;
						if (*bases == 'A' || *bases == 'a') {
							nrA++;
						} else if (*bases == 'C' || *bases == 'c') {
							nrC++;
						} else if (*bases == 'G' || *bases == 'g') {
							nrG++;
						} else if (*bases == 'T' || *bases == 't' || *bases == 'U' || *bases == 'u') {
							nrT++;
						}
					} else if (*bases == '*' || *bases == '#') {
						nrdiff++; nrdel++;
					} else if (*bases == '.' || *bases == ',') {
						sum_info[q].nummatch++;
						numbases++; nrmatch++; nummatch++;
					}
					/* other characters are ignored for counting :
						$: If this is the last position covered by the read
						> and <: Reference skip (due to CIGAR N -> RNA alignment)
					*/
					bases++;
					qs++;
				}
			}
		}
		read = DStringGetLine(&line, stdin);
		if ((int)read == -1) break;
	}
	fprintf(stdout,"#filetype\ttsv/errorrate_file\n");
	fprintf(stdout,"#fileversion\t0.99\n");
	fprintf(stdout,"#description\tinfo for estimating error rate\n");
	fprintf(stdout,"#general bam stats:\n");
	fprintf(stdout,"#nummatch\t%.0f\n",nummatch);
	fprintf(stdout,"#nummismatch\t%.0f\n",nummismatch);
	fprintf(stdout,"#numdel\t%.0f\n",numdel);
	fprintf(stdout,"#numins\t%.0f\n",numins);
	fprintf(stdout,"#mismatch_pct\t%.5f\n",100.0*nummismatch/(nummatch+nummismatch+numdel+numins));
	fprintf(stdout,"#del_pct\t%.5f\n",100.0*numdel/(nummatch+nummismatch+numdel+numins));
	fprintf(stdout,"#ins_pct\t%.5f\n",100.0*numins/(nummatch+nummismatch+numdel+numins));
	fprintf(stdout,"#error_pct\t%.5f\n",100.0*(nummismatch+numdel+numins)/(nummatch+nummismatch+numdel+numins));
	fprintf(stdout,"quality\tnummatch\tnummismatch\tnumdel\tnumins\terror_pc\n");
	for (q = 0 ; q < 60 ; q++) {
		double total;
		nummatch = sum_info[q].nummatch;
		nummismatch = sum_info[q].nummismatch;
		numdel = sum_info[q].numdel;
		numins = sum_info[q].numins;
		total = nummatch+nummismatch+numdel+numins;
		fprintf(stdout,"%d\t%.0f\t%.0f\t%.0f\t%.0f\t%.5f\t%.5f\t%.5f\t%.5f\n",
			q,nummatch,nummismatch,numdel,numins,
			total>0?100.0*nummismatch/total:0.0,
			total>0?100.0*numdel/total:0.0,
			total>0?100.0*numins/total:0.0,
			total>0?100.0*(nummismatch+numdel+numins)/total:0.0
		);
	}
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
