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

int main(int argc, char *argv[]) {
	FILE *f1,*f2,*f3;
	DStringArray *result1=NULL;
	DString *line1 = NULL;
	off_t fpos;
	uint64_t count = -1, progress = 20000000L;
	uint64_t data;
	int chr1pos,start1pos,end1pos,type1pos,max1;
	if ((argc != 8)) {
		fprintf(stderr,"Format is: bcol_indexfile file indexfile indexfilebin chrpos startpos endpos typepos");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64_or_die(argv[1],"r");
	chr1pos = atoi(argv[4]);
	start1pos = atoi(argv[5]);
	end1pos = atoi(argv[6]);
	type1pos = atoi(argv[7]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ; if (type1pos > max1) {max1 = type1pos;} ;
	f2 = fopen64_or_die(argv[2],"w");
	fprintf(f2,"# binary column\n");
	fprintf(f2,"# type wu\n");
	fprintf(f2,"chromosome\tbegin\tend\n");
	f3 = fopen64_or_die(argv[3],"w");
NODPRINT("poss: %d:%d-%d %d",chr1pos,start1pos,end1pos,type1pos)
	/* allocate */
	line1 = DStringNew();
	result1 = DStringArrayNew(max1+2);
	skip_header(f1,line1,NULL,NULL);
	fpos = ftello(f1);
	while (!DStringGetTab(line1,f1,max1,result1,1,NULL)) {
		count++;
		data = (uint64_t)fpos;
		fwrite(&data,8,1,f3);
/*
		{
		char *chromosome1;
		int start1,end1;
		chromosome1 = result1->data[chr1pos].string;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
		NODPRINT("%d\t%s\t%d\t%d\t%d",1,chromosome1,start1,end1)
		}
*/
		fpos = ftello(f1);
		if (fpos > progress) {
			fprintf(stderr,"filepos: %" PRIu64 "\n",fpos);
			progress += 20000000L;
		}
	}
	if (count == -1) {count = 0;}
	fprintf(f2,"\t0\t%" PRIu64 "\n",(count+1));
	FCLOSE(f1);
	FCLOSE(f2);
	FCLOSE(f3);
	if (line1) {DStringDestroy(line1);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
