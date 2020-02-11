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
#include <limits.h>
#include "tools.h"
#include "debug.h"


FILE *openreg(char **regions, char *prefix, char *postfix, int printheader, DString *header, DString **chromosome2, int *start2, int *end2, char *opencmd) {
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
				if (size >= 9 && strncmp(region+size-9,"unaligned",9) == 0) {
					DStringSetS(*chromosome2,"*",1);
				} else {
					DStringSetS(*chromosome2,region,size);
				}
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
	if (opencmd != NULL && opencmd[0] != '\0') {
		DStringAppend(buffer,opencmd);
		DStringAppendS(buffer," ",1);
	}
	DStringAppend(buffer,prefix);
	DStringAppendS(buffer,*regions,size);
	DStringAppend(buffer,postfix);
	if (region[size] == ' ') size++;
	*regions = region + size;
	if (opencmd == NULL || opencmd[0] == '\0') {
		o = fopen64_or_die(buffer->string,"w");
	} else {
		o = popen(buffer->string,"w");
		if (o == NULL) {
			fprintf(stderr,"error opening to command pipe: %s\n",buffer->string);
			exit(1);
		}
	}
	if (printheader && header->size) fprintf(o,"%s",header->string);
	return(o);
}

void closereg(FILE *o,char *opencmd) {
	int status;
	if (o == NULL) {return;}
	if (opencmd == NULL) {
		fclose(o);
	} else {
		status = pclose(o);
		if (status != 0)	{
			fprintf(stderr,"error closing command pipe: %s\n",opencmd);
			exit(1);
		}
	}
}

int main(int argc, char *argv[]) {
	FILE *f1,*o;
	DStringArray *result1=NULL;
	DString *header = DStringNew();
	DString *line1 = DStringNew();
	DString *chromosome1 = NULL,*chromosome2 = NULL,*curchromosome = NULL,*chromosomekeep = NULL;
	char *prefix, *postfix, *regions,*opencmd=NULL,commentchar;
	int comp,chr1pos,start1pos,end1pos,max1,printheader = 1,read,headerline;
	unsigned int numfields,pos1;
	int start1,end1,start2,end2;
	int prevstart1 = -1,prevend1 = -1,prevstart2 = -1,prevend2 = -1;
	if (argc != 10 && argc != 11) {
		fprintf(stderr,"Format is: distrreg prefix postfix printheader regions chrpos startpos endpos headerline commentchar ?opencmd?");
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
	headerline = atoi(argv[8]);
	commentchar = argv[9][0];
	if (argc == 11) opencmd = argv[10];
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	/* allocate */
	curchromosome = DStringEmtpy();
	chromosome2 = DStringNew();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(max1+2);
	/* start */
	header = DStringNew();
	read = DStringGetLine(line1, f1);
	while (read != -1) {
		if (line1->size > 0 && line1->string[0] != commentchar) break;
		DStringAppendS(header,line1->string,line1->size);
		DStringAppendS(header,"\n",1);
		read = DStringGetLine(line1, f1);
	}
	if (headerline) {
		DStringAppendS(header,line1->string,line1->size);
		DStringAppendS(header,"\n",1);
		read = DStringGetTab(line1,f1,max1,result1,0,&numfields);
	} else if (read != -1) {
		DStringSplitTab(line1,	max1, result1, 0, NULL);
	}
	o = openreg(&regions,prefix,postfix,printheader,header,&chromosome2,&start2,&end2,opencmd);
	if (o == NULL) {
		fprintf(stderr,"no regions given");
		exit(EXIT_FAILURE);
	}
	chromosomekeep = chromosome2;
	if (read != -1) do {
		pos1++;
		chromosome1 = result1->data+chr1pos;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
NODPRINT("%d\t%s\t%d\t%d",1,Loc_ChrString(chromosome1),start1,end1)
NODPRINT("%d\t%s\t%d\t%d",2,Loc_ChrString(chromosome2),start2,end2)
NODPRINT("%d\t%s\t%d\t%d",2,Loc_ChrString(curchromosome),start2,end2)
		checksortreg(curchromosome,&prevstart1,&prevend1,chromosome1,start1,end1,"stdin");
		if (o != NULL) {
			comp = DStringLocCompare(chromosome2, chromosome1);
			if (comp != 0 && chromosome2->string[chromosome2->size-1] == '_' && chromosome2->size <= chromosome1->size) {
				if (strncmp(chromosome2->string,chromosome1->string,chromosome2->size) == 0) {comp = 0;}
			}
		}
		while (o != NULL && ((comp < 0) || ((comp == 0) && (end2 <= start1)))) {
			closereg(o,opencmd);
			o = openreg(&regions,prefix,postfix,printheader,header,&chromosome2,&start2,&end2,opencmd);
			if (o == NULL) {
				chromosome2 = NULL;
				comp = -1;
				break;
			}
			comp = DStringLocCompare(chromosome2, chromosomekeep);
			if (comp != 0 && chromosome2->string[chromosome2->size-1] == '_' && chromosome2->size <= chromosomekeep->size) {
				if (strncmp(chromosome2->string,chromosomekeep->string,chromosome2->size) == 0) {comp = 0;}
			}
			if (comp < 0 || (comp == 0 && (start2 < prevstart2 || (start2 == prevstart2 && end2 < prevend2)))) {
				fprintf(stderr,"Cannot distrreg because the reglist is not correctly sorted (%s:%d-%d > %s:%d-%d)",chromosomekeep->string, prevstart2,prevend2,chromosome2->string, start2,end2);
				exit(1);
			}
			prevstart2 = start2; prevend2 = end2;
			comp = DStringLocCompare(chromosome2, chromosome1);
			if (comp != 0 && chromosome2->string[chromosome2->size-1] == '_' && chromosome2->size <= chromosome1->size) {
				if (strncmp(chromosome2->string,chromosome1->string,chromosome2->size) == 0) {comp = 0;}
			}
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
	} while (!DStringGetTab(line1,f1,max1,result1,0,&numfields));
	if (chromosomekeep != NULL) DStringDestroy(chromosomekeep);
	fclose(f1);
	while (o != NULL) {
		closereg(o,opencmd);
		o = openreg(&regions,prefix,postfix,printheader,header,&chromosome2,&start2,&end2,opencmd);
	}
	closereg(o,opencmd);
	if (line1) {DStringDestroy(line1);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
