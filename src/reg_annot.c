
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
	DString *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*data1 = NULL,*data2 = NULL,*empty=NULL;
	char *chromosome1,*chromosome2;
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,data1pos,data2pos,max1,max2;
	int nchrkeep=0,startkeep,endkeep,near,near2;
	int nchr1=0,start1,end1,nchr2=0,start2,end2;
	int error2,curchr=0,nextpos=0,datanear=-1;
	if (argc == 12) {
		datanear = atoi(argv[11]);
	} else if ((argc != 11)) {
		fprintf(stderr,"Format is: reg_annot file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2 data1pos data2pos ?datanear?");
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
	data1pos = atoi(argv[9]);
	data2pos = atoi(argv[10]);
NODPRINT("datapos %d %d",data1pos,data2pos)
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	if (data1pos > max2) {max2 = data1pos;} ; if (data2pos > max2) {max2 = data2pos;} ;
	/* allocate */
	line1 = DStringNew(); line2=DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	result1 = DStringArrayNew(max1+1);
	result2 = DStringArrayNew(max2+1);
	resultkeep = DStringArrayNew(max2+1);
	data1 = empty; data2 = empty;
	DStringGetLine(line1,f1);
	DStringGetLine(line2,f2);
	error2 = DStringGetTab(line2,f2,max2,result2);
	chromosome2 = result2[chr2pos].string;
	nchr2 = chromosomenum(chromosome2);
	sscanf(result2[start2pos].string,"%d",&start2);
	sscanf(result2[end2pos].string,"%d",&end2);
	if (data1pos != -1) {data1 = result2+data1pos;}
	if (data2pos != -1) {data2 = result2+data2pos;}
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
			nchrkeep = nchr2; startkeep = start2; endkeep = end2;
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
			if (data1pos != -1) {data1 = result2+data1pos;}
			if (data2pos != -1) {data2 = result2+data2pos;}
		}
NODPRINT("data1 %d %s",data1->size,data1->string)
NODPRINT("data2 %d %s",data2->size,data1->string)
		if (error2 || (nchr1 < nchr2) || ((nchr1 == nchr2) && (end1 <= start2))) {
			NODPRINT("no overlap")
			if (datanear != -1) {
				near = datanear;
				if (nchr1 == nchr2) {
					near = start2-end1;
					if (data1pos != -1) {data1 = result2+data1pos;}
					if (data2pos != -1) {data2 = result2+data2pos;}
				}
				if (nchr1 == nchrkeep) {
					near2 = start1-endkeep;
					if (near2 < near) {
						near = near2;
						if (data1pos != -1) {data1 = resultkeep+data1pos;}
						if (data2pos != -1) {data2 = resultkeep+data2pos;}
					}
				}
				if (near >= datanear) {
					if (data1pos == -1) {
						fprintf(stdout,"\n");
					} else if (data2pos == -1) {
						fprintf(stdout,"\t\n");
					} else {
						fprintf(stdout,"\t\t\n");
					}
				} else {
					if (data1pos == -1) {
						fprintf(stdout,"%d\n",near);
					} else if (data2pos == -1) {
						fprintf(stdout,"%s\t%d\n", data1->string,near);
					} else {
						fprintf(stdout,"%s\t%s\t%d\n", data1->string, data2->string, near);
					}
				}
			} else {
				if (data2pos == -1) {
					fprintf(stdout,"\n");
				} else {
					fprintf(stdout,"\t\n");
				}
			}
		} else {
			NODPRINT("overlap")
			if (datanear != -1) {
				near = end1-end2;
				if (near >= 0) {near = -1;}
				near2 = start2-start1-1;
				if (near2 >= 0) {near2 = -1;}
				if (near2 > near) {near = near2;}
				if (data1pos == -1) {
					fprintf(stdout,"%d\n",near);
				} else if (data2pos == -1) {
					fprintf(stdout,"%s\t%d\n", data1->string,near);
				} else {
					fprintf(stdout,"%s\t%s\t%d\n", data1->string, data2->string, near);
				}
			} else {
				if (data1pos == -1) {
					fprintf(stdout,"1\n");
				} else if (data2pos == -1) {
					fprintf(stdout,"%s\n", data1->string);
				} else {
					fprintf(stdout,"%s\t%s\n", data1->string, data2->string);
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
	exit(EXIT_SUCCESS);
}
