/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "debug.h"
#include "tools_bcol.h"

int main(int argc, char *argv[]) {
	DString *buffer = NULL, *prevchr = NULL;
	FILE *obcol;
	FILE *obin;
	char **outputa = NULL,*fieldStr,*valueStr;
	unsigned long long start;
	unsigned long long lastpos;
	DStringArray *result = NULL, *mvalues = NULL;
	DString *line = NULL,*chromosome = NULL;
	char *outfile,*type = "u",*defaultvalue = "";
	uint64_t offset, poffset = -1, size;
	int fieldcursize, valuecursize;
	int reverse = 0, isunsigned = 0, header = 0;
	int col = 0,max = 0,offsetcol = -1,chrcol = -1,mcol=0,shift,i,c,vc;
	#
	if ((argc < 4)||(argc > 11)) {
		fprintf(stderr,"Format is: bcol_make_multi output_file type mcol mvalues ?col? ?chromosomecol? ?offsetcol? ?default? ?header?\n");
		exit(EXIT_FAILURE);
	}
	outfile = argv[1];
	type = argv[2];
	mcol = atoi(argv[3]);
	max = mcol;
	mvalues = DStringArrayFromChar(argv[4],',');
	outputa = (char **)malloc(mvalues->size*sizeof(char *));
	if (argc > 5) {
		col = atoi(argv[5]);
		if (col > max) {max = col;}
	}
	if (argc > 6) {
		chrcol = atoi(argv[6]);
		if (chrcol > max) {max = chrcol;}
	}
	if (argc > 7) {
		/* this col contains the position in the chromosome */
		offsetcol = atoi(argv[7]);
		if (offsetcol > max) {max = offsetcol;}
	}
	if (argc > 8) {
		defaultvalue = argv[8];
	}
	line = DStringNew();
	if (argc > 9) {
		header = atoi(argv[9]);
		if (header) {
			skip_header(stdin,line,NULL,NULL);
		}
	}
	NODPRINT("bcol_make %s %s %d %d %d\n",outfile,type,col,chrcol,offsetcol)
	/*
		open files for writing
	 */
	start = 0;
	lastpos = -1;
	buffer = DStringNew();
	DStringAppend(buffer,outfile);
	obcol = fopen64_or_die(buffer->string,"w");
	obin = stdout;

	reverse = bcol_NeedReversing(type[0]);
	if (type[1] == 'u') {isunsigned = 1;}
	result = DStringArrayNew(max+2);
	fprintf(obcol,"# binary column\n");
	fprintf(obcol,"# type %s\n",type);
	fprintf(obcol,"# default %s\n","0");
	fprintf(obcol,"# multi %s\n",argv[4]);
	fprintf(obcol,"chromosome\tbegin\tend\n");
	poffset = -1;
	while (!DStringGetTab(line,stdin,max,result,0,NULL)) {
		if (chrcol != -1) {
			chromosome = result->data+chrcol;
			if (DStringCompare(chromosome, prevchr) != 0) {
				if (prevchr == NULL) {
					prevchr = DStringDup(chromosome);
				} else {
					if (prevchr->size > 3 && prevchr->string[0] == 'c' && prevchr->string[1] == 'h' && prevchr->string[2] == 'r') {
						shift = 3;
					} else {
						shift = 0;
					}
					fprintf(obcol,"%*.*s\t%lld\t%lld\n",prevchr->size-shift,prevchr->size-shift,prevchr->string+shift,start,start+lastpos+1);
					DStringCopy(prevchr,chromosome);
				}
				poffset = -1;
				start = 0;
				lastpos = -1;
			}
		}
		if (offsetcol != -1) {
			offset = atoll(result->data[offsetcol].string);
			if (poffset == -1) {
				start = offset;
			} else if (poffset != offset) {
				size = offset - poffset;
				if (size < 0) {
					fprintf(stderr, "error: cannot make position based bcol on unsorted data ($offset < $poffset sort on position first)\n");
					exit(EXIT_FAILURE);
				}
				while (poffset < offset) {
					for (i=0 ; i < mvalues->size ; i++) {
						bcol_printbin(obin,reverse,isunsigned,type,defaultvalue);
					}
					poffset++;
				}
				lastpos = lastpos + size;
			}
			poffset = offset+1;
		}
		NODPRINT("s=%s\n",result->data[col].string)
		for (i=0 ; i < mvalues->size ; i++) {
			outputa[i] = defaultvalue;
		}
		fieldStr = result->data[mcol].string;
		fieldcursize = 0;
		valueStr = result->data[col].string;
		valuecursize = 0;
		while (1) {
			c = fieldStr[fieldcursize];
			if (c == ',' || c =='\0' || c =='\t') {
				while(1) {
					vc = valueStr[valuecursize];
					if (vc == ',' || vc =='\0' || vc =='\t') break;
					valuecursize++;
				}
				for (i=0 ; i < mvalues->size ; i++) {
					DString *m = mvalues->data+i;
					if (fieldcursize == m->size && strncmp(m->string,fieldStr,fieldcursize) == 0) {
						outputa[i] = valueStr;
					}
				}
				if (c =='\0' || c =='\t') break;
				if (vc =='\0' || vc =='\t') break;
				fieldStr += fieldcursize + 1;
				fieldcursize = 0;
				valueStr += valuecursize + 1;
				valuecursize = 0;
			} else {
				fieldcursize++;
			}
		}
		for (i=0 ; i < mvalues->size ; i++) {
			bcol_printbin(obin,reverse,isunsigned,type,outputa[i]);
		}
		lastpos ++;
	}
	shift = 0;
	if (prevchr == NULL) {
		prevchr = DStringNew();
	} else if (prevchr->size > 3 && prevchr->string[0] == 'c' && prevchr->string[1] == 'h' && prevchr->string[2] == 'r') {
		shift = 3;
	}
	fprintf(obcol,"%*.*s\t%lld\t%lld\n",prevchr->size-shift,prevchr->size-shift,prevchr->string+shift,start,start+lastpos+1);
	if (line) {DStringDestroy(line);}
	if (result) {DStringArrayDestroy(result);}
	if (buffer) {DStringDestroy(buffer);}
	if (prevchr) {DStringDestroy(prevchr);}
	exit(EXIT_SUCCESS);
}