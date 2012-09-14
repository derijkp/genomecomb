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
	FILE *f1,*f2;
	DString *result1=NULL,*result2=NULL,*resultkeep=NULL,*resulttemp=NULL;
	DString *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*empty=NULL;
	DString **data = NULL;
	char *chromosome1,*chromosome2;
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	int datalen=0,*datapos=NULL;
	int nchrkeep=0,endkeep,near,near2;
	int nchr1=0,start1,end1,nchr2=0,start2,end2;
	int error2,curchr=0,nextpos=0,datanear=-1,i;
	if ((argc < 10)) {
		fprintf(stderr,"Format is: reg_annot file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2 datanear datapos1 ...");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	f2 = fopen64(argv[5],"r");
	chr2pos = atoi(argv[6]);
	start2pos = atoi(argv[7]);
	end2pos = atoi(argv[8]);
	datanear = atoi(argv[9]);
	datalen = argc-10;
	datapos = malloc(datalen*sizeof(int));
	data = malloc(datalen*sizeof(DString *));
	for (i = 0 ; i < datalen ; i++) {
		datapos[i] = atoi(argv[10+i]);
	}
	
NODPRINT("datapos %d %d",datapos[0],datapos[1])
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	for (i = 0 ; i < datalen ; i++) {
		if (datapos[i] > max2) {max2 = datapos[i];}
		data[i] = empty;
	}
	/* allocate */
	line1 = DStringNew(); line2=DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	result1 = DStringArrayNew(max1+1);
	result2 = DStringArrayNew(max2+1);
	resultkeep = DStringArrayNew(max2+1);
	skip_header(f1,line1);
	skip_header(f2,line2);
	error2 = DStringGetTab(line2,f2,max2,result2);
	chromosome2 = result2[chr2pos].string;
	nchr2 = chromosomenum(chromosome2);
	sscanf(result2[start2pos].string,"%d",&start2);
	sscanf(result2[end2pos].string,"%d",&end2);
	for (i = 0 ; i < datalen ; i++) {
		if (datapos[i] != -1) {data[i] = result2+datapos[i];}
	}
	while (!DStringGetTab(line1,f1,max1,result1)) {
		chromosome1 = result1[chr1pos].string;
		nchr1 = chromosomenum(chromosome1);
		sscanf(result1[start1pos].string,"%d",&start1);
		sscanf(result1[end1pos].string,"%d",&end1);
NODPRINT("%d\t%s\t%d\t%d\t%d",1,chromosome1,nchr1,start1,end1)
NODPRINT("%d\t%s\t%d\t%d\t%d",2,chromosome2,nchr2,start2,end2)
		if (nchr1 > curchr) {
			curchr = nchr1;
			nextpos = 0;
		}
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",chromosome1,start1);
			fflush(stderr);
			nextpos += 50000000;
		}
		while (!error2 && ((nchr2 < nchr1) || ((nchr2 == nchr1) && (end2 <= start1)))) {
			/* keep data of previous */
			/* to avoid allocating new memory everytime, reuse linekeep and associated data */
			nchrkeep = nchr2; endkeep = end2;
			linetemp = linekeep;
			linekeep = line2;
			line2 = linetemp;
			resulttemp = resultkeep;
			resultkeep = result2;
			result2 = resulttemp;
			/* get new line */
			error2 = DStringGetTab(line2,f2,max2,result2);
			if (error2)  {
				nchr2 = -1;
				break;
			}
			chromosome2 = result2[chr2pos].string;
			nchr2 = chromosomenum(chromosome2);
			sscanf(result2[start2pos].string,"%d",&start2);
			sscanf(result2[end2pos].string,"%d",&end2);
			for (i = 0 ; i < datalen ; i++) {
				if (datapos[i] != -1) {data[i] = result2+datapos[i];}
			}
		}
/*
DPRINT("datalen: %d",datalen)
for (i = 0; i < datalen ; i++) {
DPRINT("data[%d] %d %s",i,data[i]->size,data[i]->string)
}
*/
		if (error2 || (nchr1 < nchr2) || ((nchr1 == nchr2) && (end1 <= start2))) {
			NODPRINT("no overlap")
			if (datanear != -1) {
				near = datanear;
				if (nchr1 == nchr2) {
					near = start2-end1;
					for (i = 0 ; i < datalen ; i++) {
						if (datapos[i] != -1) {data[i] = result2+datapos[i];}
					}
				}
				if (nchr1 == nchrkeep) {
					near2 = start1-endkeep;
					if (near2 < near) {
						near = near2;
						for (i = 0 ; i < datalen ; i++) {
							if (datapos[i] != -1) {data[i] = resultkeep+datapos[i];}
						}
					}
				}
				if (near >= datanear) {
					for (i = 0 ; i < datalen ; i++) {
						fprintf(stdout,"\t");
					}
					fprintf(stdout,"\n");
				} else {
					if (datalen) {
						for (i = 0 ; i < datalen ; i++) {
							fprintf(stdout,"%s\t", data[i]->string);
						}
					}
					fprintf(stdout,"%d\n",near);
				}
			} else {
				for (i = 1 ; i < datalen ; i++) {
					fprintf(stdout,"\t");
				}
				fprintf(stdout,"\n");
			}
		} else {
			NODPRINT("overlap")
			if (datanear != -1) {
				near = end1-end2;
				if (near >= 0) {near = -1;}
				near2 = start2-start1-1;
				if (near2 >= 0) {near2 = -1;}
				if (near2 > near) {near = near2;}
				for (i = 0 ; i < datalen ; i++) {
					fprintf(stdout,"%s\t", data[i]->string);
				}
				fprintf(stdout,"%d\n",near);
			} else {
				if (!datalen) {
					fprintf(stdout,"1\n");
				} else {
					fprintf(stdout,"%s", data[0]->string);
					for (i = 1 ; i < datalen ; i++) {
						fprintf(stdout,"\t%s", data[i]->string);
					}
					fprintf(stdout,"\n");
				}
			}
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
	if (datapos) {free(datapos);}
	if (data) {free(data);}
	exit(EXIT_SUCCESS);
}
