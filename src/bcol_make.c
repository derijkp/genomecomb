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
#include <inttypes.h>
#include "tools.h"
#include "debug.h"
#include "tools_bcol.h"

int main(int argc, char *argv[]) {
	DString *buffer = NULL, *prevchr = NULL;
	FILE *obcol;
	FILE *obin;
	uint64_t start;
	DStringArray *result = NULL;
	DString *line = NULL,*chromosome = NULL;
	char *outfile,*type = "u",*defaultvalue = "",*chrname="";
	uint64_t offset, poffset, size, endpos;
	int reverse = 0, isunsigned = 0;
	int col = 0,max = 0,offsetcol = -1,endcol = -1,chrcol = -1,shift, precision = -1, i;
	if (argc != 10) {
		fprintf(stderr,"Format is: bcol_make output_file type col chromosomecol chromosomename offsetcol endcol default precision\n");
		exit(EXIT_FAILURE);
	}
	i = 1;
	outfile = argv[i++];
	type = argv[i++];
	col = atoi(argv[i++]);
	max = col;
	chrcol = atoi(argv[i++]);
	if (chrcol > max) {max = chrcol;}
	chrname = argv[i++];
	/* this col contains the position in the chromosome */
	offsetcol = atoi(argv[i++]);
	if (offsetcol > max) {max = offsetcol;}
	/* this col contains the end position in the chromosome, -1 for not used */
	endcol = atoi(argv[i++]);
	if (endcol > max) {max = endcol;}
	defaultvalue = argv[i++];
	precision = atoi(argv[i++]);
	NODPRINT(stderr,"bcol_make outfile=%s type=%s col=%d chrcol=%d chrname=%s offs=%d end=%d def=%s pr=%d\n",outfile,type,col,chrcol,chrname,offsetcol,endcol,defaultvalue,precision);
	line = DStringNew();
	/*
		open files for writing
	 */
	start = 0;
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
	if (precision != -1) {fprintf(obcol,"# precision %d\n",precision);}
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
					fprintf(obcol,"%*.*s\t%" PRId64 "\t%" PRIu64 "\n",prevchr->size-shift,prevchr->size-shift,prevchr->string+shift,start,poffset);
					fflush(obcol);
					DStringCopy(prevchr,chromosome);
				}
				poffset = -1;
				start = 0;
			}
		}
		if (offsetcol != -1) {
			offset = atoll(result->data[offsetcol].string);
			if (poffset == -1) {
				start = offset;
				poffset = offset;
			} else if (poffset != offset) {
				size = offset - poffset;
				if (size < 0) {
					fprintf(stderr, "error: cannot make position based bcol on unsorted or overlapping data ($offset < $poffset sort on position first)\n");
					exit(EXIT_FAILURE);
				}
				while (poffset < offset) {
					bcol_printbin(obin,reverse,isunsigned,type,defaultvalue);
					poffset++;
				}
			}
		} else if (poffset == -1) {
			poffset = 0;
		}
		if (endcol != -1) {
			endpos = atoll(result->data[endcol].string);
			while (poffset < endpos) {
				bcol_printbin(obin,reverse,isunsigned,type,result->data[col].string);
				poffset++;
			}
		} else {
			NODPRINT("s=%s\n",result->data[col].string)
			bcol_printbin(obin,reverse,isunsigned,type,result->data[col].string);
			poffset++;
		}
	}
	shift = 0;
	if (prevchr == NULL) {
		prevchr = DStringNewFromChar(chrname);
	} else if (prevchr->size > 3 && prevchr->string[0] == 'c' && prevchr->string[1] == 'h' && prevchr->string[2] == 'r') {
		shift = 3;
	}
	fprintf(obcol,"%*.*s\t%" PRId64 "\t%" PRIu64 "\n",prevchr->size-shift,prevchr->size-shift,prevchr->string+shift,start,poffset);
	fflush(obcol);
	if (line) {DStringDestroy(line);}
	if (result) {DStringArrayDestroy(result);}
	if (buffer) {DStringDestroy(buffer);}
	if (prevchr) {DStringDestroy(prevchr);}
	exit(EXIT_SUCCESS);
}
