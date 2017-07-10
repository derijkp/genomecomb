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
#include "gztools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	GZFILE *f2;
	FILE *f1;
	DStringArray *result1=NULL,*result2=NULL,*resultkeep=NULL,*resulttemp=NULL;
	DString *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*empty=NULL;
	DString *chromosome1 = NULL,*chromosome2 = NULL,*chromosomekeep = NULL;
	int *ontarget = NULL, *offtarget = NULL;
	int comp,chr2pos,start2pos,end2pos,max2;
	unsigned int numfields2,numfields,pos1,pos2;
	int max,depth;
	int start1,end1,start2,end2;
	int prevstart2 = -1,prevend2 = -1;
	int error2;
	if ((argc != 6)) {
		fprintf(stderr,"Format is: depth_histo file2 chrpos2 startpos2 endpos2 max");
		exit(EXIT_FAILURE);
	}
	f1 = stdin;
	if (argv[1][0] == '\0') {
		f2 = NULL;
	} else {
		f2 = gz_open(argv[1]);
	}
	chr2pos = atoi(argv[2]);
	start2pos = atoi(argv[3]);
	end2pos = atoi(argv[4]);
	max = atoi(argv[5]);
	ontarget = malloc(max*sizeof(int));
	memset(ontarget,0,max*sizeof(int));
	offtarget = malloc(max*sizeof(int));
	memset(offtarget,0,max*sizeof(int));
NODPRINT("reg_select %d %d %d %s %d %d %d %d",chr1pos,start1pos,end1pos,argv[4],chr2pos,start2pos,end2pos,datanear)
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	/* allocate */
	line1 = DStringNew(); line2=DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(2+2);
	result2 = DStringArrayNew(max2+2);
	resultkeep = DStringArrayNew(max2+2);
	if (f2 == NULL) {
		error2 = 1;
	} else {
		gz_skip_header(f2,line2,&numfields2,&pos2);
		error2 = gz_DStringGetTab(line2,f2,max2,result2,1,&numfields);	pos2++;
	}
	if (!error2) {
		check_numfieldserror(numfields,numfields2,line2,argv[4],&pos2);
		chromosome2 = result2->data+chr2pos;
		sscanf(result2->data[start2pos].string,"%d",&start2);
		sscanf(result2->data[end2pos].string,"%d",&end2);
	}
	while (!DStringGetTab(line1,f1,2,result1,0,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,3,line1,"stdin",&pos1);
		chromosome1 = result1->data+0;
		sscanf(result1->data[1].string,"%d",&end1);
		sscanf(result1->data[2].string,"%d",&depth);
		start1 = end1 - 1;
NODPRINT("%d\t%s\t%d\t%d",1,Loc_ChrString(chromosome1),start1,end1)
NODPRINT("%d\t%s\t%d\t%d",2,Loc_ChrString(chromosome2),start2,end2)
		comp = DStringLocCompare(chromosome2, chromosome1);
		while (!error2 && ((comp < 0) || ((comp == 0) && (end2 <= start1)))) {
			/* keep data of previous */
			/* to avoid allocating new memory everytime, reuse linekeep and associated data */
			chromosomekeep = chromosome2;
			linetemp = linekeep;
			linekeep = line2;
			line2 = linetemp;
			resulttemp = resultkeep;
			resultkeep = result2;
			result2 = resulttemp;
			/* get new line */
			error2 = gz_DStringGetTab(line2,f2,max2,result2,1,&numfields); pos2++;
			if (error2)  {
				chromosome2 = NULL;
				comp = -1;
				break;
			} else {
				check_numfieldserror(numfields,numfields2,line2,argv[4],&pos2);
			}
			chromosome2 = result2->data+chr2pos;
			sscanf(result2->data[start2pos].string,"%d",&start2);
			sscanf(result2->data[end2pos].string,"%d",&end2);
			
			comp = DStringLocCompare(chromosome2, chromosomekeep);
			if (comp < 0 || (comp == 0 && (start2 < prevstart2 || (start2 == prevstart2 && end2 < prevend2)))) {
				fprintf(stderr,"Cannot annotate because the database file is not correctly sorted (sort correctly using \"cg select -s -\")");
				exit(1);
			}
			prevstart2 = start2; prevend2 = end2;
			comp = DStringLocCompare(chromosome2, chromosome1);
		}
		if (depth > max) {depth=max;}
		if (error2 || (comp > 0) || ((comp == 0) && (end1 <= start2 && start1 != start2))) {
			NODPRINT("no overlap")
			offtarget[depth]++;
		} else {
			NODPRINT("overlap")
			ontarget[depth]++;
		}
	}
	fclose(f1);
	if (f2 != NULL) {
		gz_close(f2);
	}
	fprintf(stdout,"depth\tontarget\tofftarget\n");
	for(depth = 1; depth <= max; depth++) {
		if (ontarget[depth] > 0 || offtarget[depth] > 0) fprintf(stdout,"%d\t%d\t%d\n",depth,ontarget[depth],offtarget[depth]);
	}
	free(ontarget); free(offtarget);
	if (line1) {DStringDestroy(line1);}
	if (line2) {DStringDestroy(line2);}
	if (linekeep) {DStringDestroy(linekeep);}
	if (empty) {DStringDestroy(empty);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	if (resultkeep) {DStringArrayDestroy(resultkeep);}
	exit(EXIT_SUCCESS);
}
