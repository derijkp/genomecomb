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
	FILE *f1,*f2;
	DStringArray *result1=NULL,*result2=NULL,*resultkeep=NULL,*resulttemp=NULL;
	DString *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*empty=NULL;
	DString *chromosome1 = NULL,*chromosome2 = NULL,*curchromosome = NULL,*chromosomekeep = NULL;
	int comp,chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	unsigned int numfields1,numfields2,numfields,pos1,pos2;
	int endkeep=-1,near,near2;
	int start1,end1,start2,end2;
	int prevstart1 = -1,prevend1 = -1,prevstart2 = -1,prevend2 = -1;
	int error2,datanear=-1;
	if ((argc != 9)) {
		fprintf(stderr,"Format is: reg_select chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2 datanear");
		exit(EXIT_FAILURE);
	}
	f1 = stdin;
	chr1pos = atoi(argv[1]);
	start1pos = atoi(argv[2]);
	end1pos = atoi(argv[3]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	f2 = fopen64_or_die(argv[4],"r");
	chr2pos = atoi(argv[5]);
	start2pos = atoi(argv[6]);
	end2pos = atoi(argv[7]);
	datanear = atoi(argv[8]);
NODPRINT("reg_select %d %d %d %s %d %d %d %d",chr1pos,start1pos,end1pos,argv[4],chr2pos,start2pos,end2pos,datanear)
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	/* allocate */
	line1 = DStringNew(); line2=DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	curchromosome = DStringEmtpy();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(max1+2);
	result2 = DStringArrayNew(max2+2);
	resultkeep = DStringArrayNew(max2+2);
	skip_header(f1,line1,&numfields1,&pos1);
	fprintf(stdout,"%s\n",line1->string);
	skip_header(f2,line2,&numfields2,&pos2);
	error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields);	pos2++;
	if (!error2) {
		check_numfieldserror(numfields,numfields2,line2,argv[4],&pos2);
		chromosome2 = result2->data+chr2pos;
		sscanf(result2->data[start2pos].string,"%d",&start2);
		sscanf(result2->data[end2pos].string,"%d",&end2);
	}
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
		while (!error2 && ((comp < 0) || ((comp == 0) && (end2 <= start1)))) {
			/* keep data of previous */
			/* to avoid allocating new memory everytime, reuse linekeep and associated data */
			chromosomekeep = chromosome2; endkeep = end2;
			linetemp = linekeep;
			linekeep = line2;
			line2 = linetemp;
			resulttemp = resultkeep;
			resultkeep = result2;
			result2 = resulttemp;
			/* get new line */
			error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields); pos2++;
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
		if (error2 || (comp > 0) || ((comp == 0) && (end1 <= start2 && start1 != start2))) {
			NODPRINT("no overlap")
			if (datanear != -1) {
				near = datanear;
				if (comp == 0) {
					near = start2-end1;
				}
				if (DStringLocCompare(chromosome1,chromosomekeep) == 0) {
					near2 = start1-endkeep;
					if (near2 < near) {
						near = near2;
					}
				}
				if (near < datanear) {
					fprintf(stdout,"%s\n",line1->string);
				}
			}
		} else {
			NODPRINT("overlap")
			fprintf(stdout,"%s\n",line1->string);
		}
	}
	fclose(f1);
	fclose(f2);
	if (line1) {DStringDestroy(line1);}
	if (line2) {DStringDestroy(line2);}
	if (linekeep) {DStringDestroy(linekeep);}
	if (empty) {DStringDestroy(empty);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	if (resultkeep) {DStringArrayDestroy(resultkeep);}
	exit(EXIT_SUCCESS);
}
