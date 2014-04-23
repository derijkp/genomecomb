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
	FILE *f1=stdin,*f2;
	DStringArray *result1=NULL,*result2=NULL,*resultkeep=NULL,*resulttemp=NULL;
	DString *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*empty=NULL,*cigar=NULL,*qual=NULL,*seq=NULL;
	DString *chromosome1 = NULL,*chromosome2 = NULL,*curchromosome = NULL,*chromosomekeep = NULL;
	char *cur;
	ssize_t read;
	unsigned int curpos=0;
	int comp,chr2pos,start2pos,end2pos,outerstart2pos,outerend2pos,max2;
	unsigned int numfields2,numfields,pos1,pos2;
	int endkeep=-1,count;
	int start1,end1,start2,end2,outerstart2,outerend2;
	int prevstart1 = -1, prevend1 = -1, prevstart2 = -1,prevend2 = -1;
	int error2;
	if ((argc != 7)) {
		fprintf(stderr,"Format is: sam_clipamplicons ampliconsfile chrpos startpos endpos outerstartpos outerendpos");
		exit(EXIT_FAILURE);
	}
	f2 = fopen64_or_die(argv[1],"r");
	chr2pos = atoi(argv[2]);
	start2pos = atoi(argv[3]);
	end2pos = atoi(argv[4]);
	outerstart2pos = atoi(argv[5]);
	outerend2pos = atoi(argv[6]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	if (outerstart2pos > max2) {max2 = outerstart2pos;} ; if (outerend2pos > max2) {max2 = outerend2pos;} ;
	/* allocate */
	line1 = DStringNew(); line2=DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	curchromosome = DStringEmtpy();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(10+2);
	result2 = DStringArrayNew(max2+2);
	resultkeep = DStringArrayNew(max2+2);
	read = DStringGetLine(line1, f1);
	while (read != -1) {
		if (line1->string[0] != '@') break;
		fprintf(stdout,"%s\n",line1->string);
		read = DStringGetLine(line1, f1); curpos++;
		pos1++;
	}
	skip_header(f2,line2,&numfields2,&pos2);
	error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields);	pos2++;
	if (!error2) {
		check_numfieldserror(numfields,numfields2,line2,argv[1],&pos2);
		chromosome2 = result2->data+chr2pos;
		sscanf(result2->data[start2pos].string,"%d",&start2);
		sscanf(result2->data[end2pos].string,"%d",&end2);
		sscanf(result2->data[start2pos].string,"%d",&outerstart2);
		sscanf(result2->data[outerend2pos].string,"%d",&outerend2);
	}
	DStringSplitTab(line1,10,result1,0,&numfields);
	while (1) {
		pos1++;
		/* check_numfieldserror(numfields,15,line1,"stdin",&pos1); */
		chromosome1 = result1->data+2;
		sscanf(result1->data[3].string,"%d",&start1);
		end1 = start1;
		cigar = result1->data+5;
		seq = result1->data+9;
		qual = result1->data+10;
		if (chromosome1->string[0] == '*') {
			NODPRINT("unmapped")
			fprintf(stdout,"%s\n",line1->string); 
		} else {
			checksortreg(curchromosome,&prevstart1,&prevend1,chromosome1,start1,end1,"stdin");
			comp = DStringLocCompare(chromosome2, chromosome1);
			while (!error2 && ((comp < 0) || ((comp == 0) && (start1 < outerstart2)))) {
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
					check_numfieldserror(numfields,numfields2,line2,argv[1],&pos2);
				}
				chromosome2 = result2->data+chr2pos;
				sscanf(result2->data[start2pos].string,"%d",&start2);
				sscanf(result2->data[end2pos].string,"%d",&end2);
				sscanf(result2->data[start2pos].string,"%d",&outerstart2);
				sscanf(result2->data[outerend2pos].string,"%d",&outerend2);
				
				comp = DStringLocCompare(chromosome2, chromosomekeep);
				if (comp < 0 || (comp == 0 && (start2 < prevstart2 || (start2 == prevstart2 && end2 < prevend2)))) {
					fprintf(stderr,"Cannot annotate because the database file is not correctly sorted (sort correctly using \"cg select -s -\")");
					exit(1);
				}
				prevstart2 = start2; prevend2 = end2;
				comp = DStringLocCompare(chromosome2, chromosome1);
			}
			if (comp == 0) {
				NODPRINT("overlap")
			}
			cur = seq->string;
			if (*cur != '*') {
				count = 30;
				while (count--) {*cur++ = 'N';}
			}
			cur = qual->string;
			if (*cur != '*') {
				count = 30;
				while (count--) {*cur++ = '!';}
			}
			fprintf(stdout,"%s\n",line1->string);
/*
			if (error2 || (comp > 0) || ((comp == 0) && (end1 <= start2))) {
				NODPRINT("no overlap")
				fprintf(stdout,"%s\n",line1->string);
			} else {
				int i
				NODPRINT("overlap")
				fprintf(stdout,"%s\n",line1->string);
				fprintf(stdout,"%*.*s\t%*.*s\t%*.*s\t%d\t%*.*s\t%*.*s\t%*.*s\n",
					result1->data[0].size,result1->data[0].size,result1->data[0].string,
					result1->data[1].size,result1->data[1].size,result1->data[1].string,
					chromosome1->size,chromosome1->size,chromosome1->string,
					start1,
					result1->data[4].size,result1->data[4].size,result1->data[4].string,
					cigar->size,cigar->size,cigar->string,
					result1->data[6].size,result1->data[6].size,result1->data[6].string
				);
			}
*/
		}
		if (DStringGetTab(line1,f1,10,result1,0,&numfields)) break;
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
