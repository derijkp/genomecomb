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
	DStringArray *result1=NULL,*result2=NULL,*resultkeep=NULL,*resulttemp=NULL;
	DString *line1 = NULL,*line2 = NULL,*linekeep = NULL,*linetemp = NULL,*empty=NULL;
	DString **data = NULL;
	DString *chromosome1 = NULL,*chromosome2 = NULL,*curchromosome = NULL,*chromosomekeep = NULL;
	DString *prevchromosome1 = NULL, *prevchromosome2 = NULL;
	int comp,chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	int datalen=0,*datapos=NULL;
	unsigned int numfields1,numfields2,numfields,pos1=0,pos2=0;
	int endkeep=-1,near,near2;
	int start1,end1,start2,end2;
	int prevstart1 = -1,prevend1 = -1,prevstart2 = -1,prevend2 = -1;
	int error2,nextpos=0,datanear=-1,i;
	if ((argc < 10)) {
		fprintf(stderr,"Format is: reg_annot file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2 datanear datapos1 ...");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64_or_die(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	if (strlen(argv[5]) == 1 && argv[5][0] == '-') {
		f2 = stdin;
	} else {
		f2 = fopen64_or_die(argv[5],"r");
	}
	chr2pos = atoi(argv[6]);
	start2pos = atoi(argv[7]);
	end2pos = atoi(argv[8]);
	datanear = atoi(argv[9]);
NODPRINT("reg_annot %s %d %d %d %s %d %d %d %d",argv[1],chr1pos,start1pos,end1pos,argv[5],chr2pos,start2pos,end2pos,datanear)
	datalen = argc-10;
	datapos = malloc(datalen*sizeof(int));
	data = malloc(datalen*sizeof(DString *));
	for (i = 0 ; i < datalen ; i++) {
		datapos[i] = atoi(argv[10+i]);
NODPRINT("%d",datapos[i])
	}
	
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	for (i = 0 ; i < datalen ; i++) {
		if (datapos[i] > max2) {max2 = datapos[i];}
		data[i] = empty;
	}
	/* allocate */
	line1 = DStringNew(); line2=DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	prevchromosome1 = DStringNew();	prevchromosome2 = DStringNew();
	curchromosome = DStringEmtpy();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(max1+2);
	result2 = DStringArrayNew(max2+2);
	resultkeep = DStringArrayNew(max2+2);
	skip_header(f1,line1,&numfields1,&pos1);
	skip_header(f2,line2,&numfields2,&pos2);
	error2 = DStringGetTab(line2,f2,max2,result2,1,&numfields);	pos2++;
	if (!error2) {
		check_numfieldserror(numfields,numfields2,line2,argv[5],&pos2);
		chromosome2 = result2->data+chr2pos;
		sscanf(result2->data[start2pos].string,"%d",&start2);
		sscanf(result2->data[end2pos].string,"%d",&end2);
	}
	for (i = 0 ; i < datalen ; i++) {
		if (datapos[i] != -1) {data[i] = result2->data+datapos[i];}
	}
	while (!DStringGetTab(line1,f1,max1,result1,1,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,numfields1,line1,argv[1],&pos1);
		chromosome1 = result1->data+chr1pos;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
NODPRINT("%d\t%s\t%d\t%d",1,Loc_ChrString(chromosome1),start1,end1)
NODPRINT("%d\t%s\t%d\t%d",2,Loc_ChrString(chromosome2),start2,end2)
NODPRINT("%d\t%s\t%d\t%d",2,Loc_ChrString(curchromosome),start2,end2)
		if (checksortreg(curchromosome,&prevstart1,&prevend1,chromosome1,start1,end1,argv[1])) {
			nextpos = 0;
		}
		/* if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",Loc_ChrString(chromosome1),start1);
			fflush(stderr);
			nextpos += 50000000;
		} */
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
				check_numfieldserror(numfields,numfields2,line2,argv[5],&pos2);
			}
			chromosome2 = result2->data+chr2pos;
			sscanf(result2->data[start2pos].string,"%d",&start2);
			sscanf(result2->data[end2pos].string,"%d",&end2);
			checksortreg(prevchromosome2,&prevstart2,&prevend2,chromosome2,start2,end2,"database file");
			for (i = 0 ; i < datalen ; i++) {
				if (datapos[i] != -1) {data[i] = result2->data+datapos[i];}
			}
			comp = DStringLocCompare(chromosome2, chromosome1);
		}
/*
DPRINT("datalen: %d",datalen)
for (i = 0; i < datalen ; i++) {
DPRINT("data[%d] %d %s",i,data[i]->size,data[i]->string)
}
*/
		if (error2 || (comp > 0) || ((comp == 0) && (end1 <= start2))) {
			NODPRINT("no overlap")
			if (datanear != -1) {
				near = datanear;
				if (comp == 0) {
					near = start2-end1;
					for (i = 0 ; i < datalen ; i++) {
						if (datapos[i] != -1) {data[i] = result2->data+datapos[i];}
					}
				}
				if (DStringLocCompare(chromosome1,chromosomekeep) == 0) {
					near2 = start1-endkeep;
					if (near2 < near) {
						near = near2;
						for (i = 0 ; i < datalen ; i++) {
							if (datapos[i] != -1) {data[i] = resultkeep->data+datapos[i];}
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
							if (data[i] == NULL) {
								fprintf(stdout,"\t");
							} else {
								fprintf(stdout,"%s\t", data[i]->string);
							}
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
					if (data[i] == NULL) {
						fprintf(stdout,"\t");
					} else {
						fprintf(stdout,"%s\t", data[i]->string);
					}
				}
				fprintf(stdout,"%d\n",near);
			} else {
				if (!datalen) {
					fprintf(stdout,"1\n");
				} else {
					fprintf(stdout,"%s", data[0]->string);
					for (i = 1 ; i < datalen ; i++) {
						if (data[i] == NULL) {
							fprintf(stdout,"\t");
						} else {
							fprintf(stdout,"\t%s", data[i]->string);
						}
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
	if (curchromosome) {DStringDestroy(curchromosome);}
	if (prevchromosome1) {DStringDestroy(prevchromosome1);}
	if (prevchromosome2) {DStringDestroy(prevchromosome2);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	if (resultkeep) {DStringArrayDestroy(resultkeep);}
	if (datapos) {free(datapos);}
	if (data) {free(data);}
	exit(EXIT_SUCCESS);
}
