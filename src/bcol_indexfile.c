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

int main(int argc, char *argv[]) {
	FILE *f1,*f2,*f3;
	DStringArray *result1=NULL;
	DString *line1 = NULL;
	off_t fpos;
	uint64_t count = -1,next = 4294967296LL, offset = 0L, progress = 50000000L;
	uint32_t data;
	int chr1pos,start1pos,end1pos,type1pos,max1;
	if ((argc != 8)) {
		fprintf(stderr,"Format is: bcol_indexfile file indexfile indexfilebin chrpos startpos endpos typepos");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[4]);
	start1pos = atoi(argv[5]);
	end1pos = atoi(argv[6]);
	type1pos = atoi(argv[7]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ; if (type1pos > max1) {max1 = type1pos;} ;
	f2 = fopen64(argv[2],"w");
	fprintf(f2,"# binary column\n");
	fprintf(f2,"# type lineindex\n");
	fprintf(f2,"begin\ttype\toffset\n");
	f3 = fopen64(argv[3],"w");
NODPRINT("poss: %d:%d-%d %d",chr1pos,start1pos,end1pos,type1pos)
	/* allocate */
	line1 = DStringNew();
	result1 = DStringArrayNew(max1+2);
	skip_header(f1,line1);
	fpos = ftello(f1);
	if (fpos >= next) {
		fprintf(stderr,"Header too long");
		exit(EXIT_FAILURE);
	}
	fprintf(f2,"%d\t%s\t%d\n",0,"iu",0);
	while (!DStringGetTab(line1,f1,max1,result1,1)) {
		count++;
		data = (uint32_t)(fpos-offset);
		fwrite(&data,4,1,f3);
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
			fprintf(stderr,"filepos: %llu\n",(unsigned long long)fpos);
			progress += 50000000L;
		}
		if (fpos >= next) {
			while (fpos >= next) {
				offset=next;
				next = next << 1;
			}
			fprintf(f2,"%ju\t%s\t%ju\n",(uintmax_t)count,"iu",(uintmax_t)offset);
		}
	}
	fprintf(f2,"%ju\t%s\t%ju\n",(uintmax_t)count,"end",(uintmax_t)offset);
	fclose(f1);
	fclose(f2);
	fclose(f3);
	if (line1) {DStringDestroy(line1);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
