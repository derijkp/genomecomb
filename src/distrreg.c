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
#include "tools.h"
#include "gztools.h"
#include "debug.h"


FILE *openreg(char **regions, char *prefix, char *postfix, int printheader, DString *header, DString **chromosome2, int *start2, int *end2) {
	FILE *o = NULL;
	DString *buffer = DStringNew();
	char *region = *regions;
	int size = 0, step = 1;
	if (*region == '\0') return(NULL);
	*chromosome2 = DStringNew();
	*start2 = 0 ; *end2 = INT_MAX;
	while(1) {
		if (region[size] == '\0' || region[size] == ' ' || region[size] == ':' || region[size] == '-') {
			if (step == 1) {
				DStringSetS(*chromosome2,region,size);
				if (region[size] == '\0' || region[size] == ' ') break;
				*start2 = atoi(region+size+1);
				step = 2;
			} else if (step == 2) {
				if (region[size] == '\0' || region[size] == ' ') break;
				*end2 = atoi(region+size+1);
				if (*end2 < *start2) {
					fprintf(stderr,"error parsing region: end %d < begin %d\n",*end2,*start2);
					exit(EXIT_FAILURE);
				}
				step = 3;
			} else if (region[size] == '\0' || region[size] == ' ') {
				break;
			} else {
				fprintf(stderr,"error parsing region \"%30.30s...\" at \"%30.30s...\"\n",region,region+size);
				exit(EXIT_FAILURE);
			}
			
		}
		size++;
	}
	DStringAppend(buffer,prefix);
	DStringAppendS(buffer,*regions,size);
	DStringAppend(buffer,postfix);
	if (region[size] == ' ') size++;
	*regions = region + size;
	o = fopen64_or_die(buffer->string,"a");
	if (printheader) fprintf(o,"%s\n",header->string);
	return(o);
}

int main(int argc, char *argv[]) {
	FILE *f1,*o;
	DStringArray *result1=NULL;
	DString *header;
	DString *line1 = NULL;
	DString *chromosome1 = NULL,*chromosome2 = NULL,*curchromosome = NULL,*chromosomekeep = NULL;
	char *prefix, *postfix, *regions;
	int comp,chr1pos,start1pos,end1pos,max1,printheader = 1;
	unsigned int numfields1,numfields,pos1;
	int start1,end1,start2,end2;
	int prevstart1 = -1,prevend1 = -1,prevstart2 = -1,prevend2 = -1;
	if ((argc != 8)) {
		fprintf(stderr,"Format is: distrreg prefix postfix printheader regions chrpos startpos endpos");
		exit(EXIT_FAILURE);
	}
	f1 = stdin;
	prefix = argv[1];
	postfix = argv[2];
	printheader = atoi(argv[3]);
	regions = argv[4];
	chr1pos = atoi(argv[5]);
	start1pos = atoi(argv[6]);
	end1pos = atoi(argv[7]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	/* allocate */
	line1 = DStringNew();
	curchromosome = DStringEmtpy();
	chromosome2 = DStringNew();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(max1+2);
	/* start */
	header = DStringNew();
	skip_header(f1,header,&numfields1,&pos1);
	o = openreg(&regions,prefix,postfix,printheader,header,&chromosome2,&start2,&end2);
	if (o == NULL) {
		fprintf(stderr,"no regions given");
		exit(EXIT_FAILURE);
	}
	chromosomekeep = chromosome2;
	while (!DStringGetTab(line1,f1,max1,result1,0,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,numfields1,line1,"stdin",&pos1);
		chromosome1 = result1->data+chr1pos;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
NODPRINT("%d\t%s\t%d\t%d",1,Loc_ChrString(chromosome1),start1,end1)
NODPRINT("%d\t%s\t%d\t%d",2,Loc_ChrString(chromosome2),start2,end2)
NODPRINT("%d\t%s\t%d\t%d",2,Loc_ChrString(curchromosome),start2,end2)
		checksortreg(curchromosome,&prevstart1,&prevend1,chromosome1,start1,end1,"stdin");
		comp = DStringLocCompare(chromosome2, chromosome1);
		while (o != NULL && ((comp < 0) || ((comp == 0) && (end2 <= start1)))) {
			fclose(o);
			o = openreg(&regions,prefix,postfix,printheader,header,&chromosome2,&start2,&end2);
			if (o == NULL) {
				chromosome2 = NULL;
				comp = -1;
				break;
			}
			comp = DStringLocCompare(chromosome2, chromosomekeep);
			if (comp < 0 || (comp == 0 && (start2 < prevstart2 || (start2 == prevstart2 && end2 < prevend2)))) {
				fprintf(stderr,"Cannot distrreg because the reglist is not correctly sorted (%s:%d-%d > %s:%d-%d)",chromosomekeep->string, prevstart2,prevend2,chromosome2->string, start2,end2);
				exit(1);
			}
			prevstart2 = start2; prevend2 = end2;
			comp = DStringLocCompare(chromosome2, chromosome1);
			if (chromosomekeep != NULL) DStringDestroy(chromosomekeep);
			chromosomekeep = chromosome2;
		}
		if (o == NULL || (comp > 0) || ((comp == 0) && (end1 <= start2 && start1 != start2))) {
			fprintf(stderr,"variants outside of distrreg regions:\n%s\n",line1->string);
			/* exit(EXIT_FAILURE); */
		} else {
			NODPRINT("overlap")
			fprintf(o,"%s\n",line1->string);
		}
	}
	if (chromosomekeep != NULL) DStringDestroy(chromosomekeep);
	fclose(f1);
	while (o != NULL) {
		fclose(o);
		o = openreg(&regions,prefix,postfix,printheader,header,&chromosome2,&start2,&end2);
	}
	if (line1) {DStringDestroy(line1);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}

