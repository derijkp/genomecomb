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
	FILE *f1;
	DStringArray *result1=NULL;
	DString *line1 = NULL;
	DString *prevchromosome1 = DStringNew();
	DString *prevtype1 = DStringNew();
	DString *prevalt1 = DStringNew();
	DString *chromosome1=NULL,*type1 = NULL,*alt1 = NULL;
	unsigned int numfields1,numfields,pos1;
	int prevstart1 = -1,prevend1 = -1;
	int chr1pos,start1pos,end1pos,type1pos,alt1pos,max1;
	int comp,comptype,compalt;
	int start1,end1;
	if ((argc != 6)) {
		fprintf(stderr,"Format is: check_sort chrpos1 startpos1 endpos1 type1pos alt1pos");
		exit(EXIT_FAILURE);
	}
	f1 = stdin;
	chr1pos = atoi(argv[1]);
	start1pos = atoi(argv[2]);
	end1pos = atoi(argv[3]);
	type1pos = atoi(argv[4]);
	alt1pos = atoi(argv[5]);
	max1 = chr1pos;
	if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	if (type1pos > max1) {max1 = type1pos;} ; if (alt1pos > max1) {max1 = alt1pos;} ;
	line1 = DStringNew();
	/* The following allocation is not destroyed at end as it may point to something else */
	/* This will leak mem, but as the prog is finished anyway ... */
	result1 = DStringArrayNew(max1+2);
	skip_header(f1,line1,&numfields1,&pos1);
	while (!DStringGetTab(line1,f1,max1,result1,1,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,numfields1,line1,argv[1],&pos1);
		chromosome1 = result1->data+chr1pos;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
		type1 = result1->data+type1pos;
		alt1 = result1->data+alt1pos;
NODPRINT("line1 (a=%3.3s) %s,%d,%d %s",type1->string,chromosome1->string,start1,end1,alt1->string)
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,Loc_ChrString(chromosome1),start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,Loc_ChrString(chromosome2),start2,end2);
*/
	 	comp = DStringLocCompare(chromosome1,prevchromosome1);
		if (type1pos != -1) {
			comptype = DStringLocCompare(type1,prevtype1);
		} else {
			comptype = 0;
		}
		if (alt1pos != -1) {
			compalt = DStringLocCompare(alt1,prevalt1);
		} else {
			compalt = 0;
		}
		if (comp < 0 || (comp == 0 && 
			(start1 < prevstart1 || (start1 == prevstart1 && 
			(end1 < prevend1 || (end1 == prevend1 &&
			(comptype < 0 || (comptype == 0 && compalt < 0)
		))))))) {
			fprintf(stderr,"file is not correctly sorted (sort correctly using \"cg select -s -\")\n");
			if (alt1pos != -1) {
				fprintf(stderr,"%s:%d-%d:%s:%s came before %s:%d-%d:%s:%s\n",prevchromosome1->string,prevstart1,prevend1,prevtype1->string,prevalt1->string, chromosome1->string,start1,end1,type1->string,alt1->string);
			} else if (type1pos != -1) {
				fprintf(stderr,"%s:%d-%d:%s came before %s:%d-%d:%s\n",prevchromosome1->string,prevstart1,prevend1,prevtype1->string, chromosome1->string,start1,end1,type1->string);
			} else {
				fprintf(stderr,"%s:%d-%d came before %s:%d-%d\n",prevchromosome1->string,prevstart1,prevend1, chromosome1->string,start1,end1);
			}
			exit(1);
		} else if (comp > 0) {
			DStringCopy(prevchromosome1,chromosome1);
		}
		if (alt1pos != -1 && compalt != 0) {
			DStringCopy(prevalt1,alt1);
		}
		if (type1pos != -1 && comptype != 0) {
			DStringCopy(prevtype1,type1);
		}
		prevstart1 = start1; prevend1 = end1;
	}
	fclose(f1);
	DStringDestroy(prevchromosome1); DStringDestroy(prevtype1); DStringDestroy(prevalt1);
	if (line1) {DStringDestroy(line1);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
