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
#include <errno.h>
#include "tools.h"
#include "debug.h"

int multireg_next(FILE *f1,DString *line1,int chr1pos, int start1pos, int end1pos, int max1,DString *result1,int *nchr1,int *start1,int *end1,DString **cur1) {
	int error1;
	error1 = DStringGetTab(line1,f1,max1,result1);
	if (error1) {goto error;}
	*nchr1 = chromosomenum(result1[chr1pos].string);
	error1 = sscanf(result1[start1pos].string,"%d",start1);
	if (!error1) {goto error;}
	error1 = sscanf(result1[end1pos].string,"%d",end1);
	if (!error1) {goto error;}
	if (cur1 != NULL) {
		*cur1 = result1+max1+1;
	}
	return 0;
	error:
		*nchr1 = NOCHROM;
		*start1 = -1;
		*end1 = -1;
		return 1;
}

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	DString *result1=NULL,*result2=NULL;
	DString *line1 = NULL,*line2 = NULL,*cur1 = NULL;
	char chromosome[10],*nulldata;
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	int nchr1=0,start1=-1,end1=-1,nchr2=0,start2=-1,end2=-1;
	int error1,error2,curchr=0,nextpos=0;
	if ((argc != 10)) {
		fprintf(stderr,"Format is: reg_annot file1 chrpos1 startpos1 endpos1 nulldata file2 chrpos2 startpos2 endpos2\n");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	if (chr1pos != 0 || start1pos != 1 || end1pos != 2) {
		fprintf(stderr,"file1 must (currently) contain chrpos1,startpos1,endpos1 as first columns\n");
		exit(EXIT_FAILURE);
	}
	nulldata = argv[5];
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	f2 = fopen64(argv[6],"r");
	chr2pos = atoi(argv[7]);
	start2pos = atoi(argv[8]);
	end2pos = atoi(argv[9]);
NODPRINT("%d %d %d ""%s"" %d %d %d",chr1pos,start1pos,end1pos,nulldata,chr2pos,start2pos,end2pos)
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	/* allocate */
	line1 = DStringNew(); line2=DStringNew();
	result1 = DStringArrayNew(max1+2);
	result2 = DStringArrayNew(max2+2);
	skip_header(f1,line1);
	skip_header(f2,line2);
	error1 = multireg_next(f1,line1,chr1pos,start1pos,end1pos,max1,result1,&nchr1,&start1,&end1,&cur1);
	error2 = multireg_next(f2,line2,chr2pos,start2pos,end2pos,max2,result2,&nchr2,&start2,&end2,NULL);
	while (!error1 || !error2) {
NODPRINT("%d:%d-%d %s",nchr1,start1,end1,cur1->string)
NODPRINT("%d:%d-%d",nchr2,start2,end2)
		if (nchr1 > curchr && nchr2 > curchr) {curchr = (nchr1<nchr2)?nchr1:nchr2; nextpos = 0;}
		if (start1 >= nextpos || start2 >= nextpos) {
			fprintf(stderr, "%s-%d\n",num2chromosome(chromosome,(nchr1<nchr2)?nchr1:nchr2),nextpos);
			fflush(stderr);
			nextpos += 10000000;
		}
		if ((nchr1 == NOCHROM) && (nchr2 == NOCHROM)) break;
		if ((nchr2 < nchr1) || (nchr1 == nchr2 && end2 <= start1)) {
			fprintf(stdout,"%s\t%d\t%d\t%s\t%d\n",num2chromosome(chromosome,nchr2),start2,end2,nulldata,1);
			error2 = multireg_next(f2,line2,chr2pos,start2pos,end2pos,max2,result2,&nchr2,&start2,&end2,NULL);
		} else if ((nchr1 < nchr2) || (nchr1 == nchr2 && end1 <= start2)) {
			fprintf(stdout,"%s\t%d\t%d\t%s\t%d\n",num2chromosome(chromosome,nchr1),start1,end1,cur1->string,0);
			error1 = multireg_next(f1,line1,chr1pos,start1pos,end1pos,max1,result1,&nchr1,&start1,&end1,&cur1);
		} else {
			if (start1 < start2) {
				fprintf(stdout,"%s\t%d\t%d\t%s\t%d\n",num2chromosome(chromosome,nchr1),start1,start2,cur1->string,0);
				start1 = start2;
			}
			if (start2 < start1) {
				fprintf(stdout,"%s\t%d\t%d\t%s\t%d\n",num2chromosome(chromosome,nchr1),start2,start1,nulldata,1);
				start2 = start1;
			}
			if (end2 < end1) {
				fprintf(stdout,"%s\t%d\t%d\t%s\t%d\n",num2chromosome(chromosome,nchr1),start1,end2,cur1->string,1);
				start1 = end2;
				error2 = multireg_next(f2,line2,chr2pos,start2pos,end2pos,max2,result2,&nchr2,&start2,&end2,NULL);
			} else if (end1 < end2) {
				fprintf(stdout,"%s\t%d\t%d\t%s\t%d\n",num2chromosome(chromosome,nchr1),start1,end1,cur1->string,1);
				start2 = end1;
				error1 = multireg_next(f1,line1,chr1pos,start1pos,end1pos,max1,result1,&nchr1,&start1,&end1,&cur1);
			} else {
				fprintf(stdout,"%s\t%d\t%d\t%s\t%d\n",num2chromosome(chromosome,nchr1),start1,end1,cur1->string,1);
				error1 = multireg_next(f1,line1,chr1pos,start1pos,end1pos,max1,result1,&nchr1,&start1,&end1,&cur1);
				error2 = multireg_next(f2,line2,chr2pos,start2pos,end2pos,max2,result2,&nchr2,&start2,&end2,NULL);
			}
		}
	}
	fclose(f1);
	fclose(f2);
	if (line1) {DStringDestroy(line1);}
	if (line2) {DStringDestroy(line2);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	exit(EXIT_SUCCESS);
}
