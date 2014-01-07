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
#include "tools_bcol.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	FILE *f1;
	BCol *fbcol = NULL;
	DStringArray *result1=NULL;
	DString *line1 = NULL,*linekeep = NULL,*empty=NULL;
	DString *chromosome1 = NULL,*chromosome2 = NULL,*curchromosome = NULL,*chromosomekeep = NULL;
	char *bcolfile;
	unsigned int numfields1,numfields,pos1;
	int comp,chr1pos,start1pos,end1pos,max1,loaded2 = 0;
	int bcolpos;
	int start1,end1;
	int prevstart1 = -1,prevend1 = -1;
	int nextpos=0,c;
	if ((argc < 6)) {
		fprintf(stderr,"Format is: bcol_annot file1 chrpos1 startpos1 endpos1 bcol1chr bcol1file ?bcol2chr bcol2file? ...\n");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64_or_die(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	bcolpos = 5; argc -= 4;
	chromosome2 = DStringNew();
	DStringSet(chromosome2,argv[bcolpos]);
	bcolfile = argv[bcolpos+1];
	NODPRINT("reg_annot %s %d %d %d\n",argv[1],chr1pos,start1pos,end1pos)
	/* allocate */
	line1 = DStringNew(); linekeep=DStringNew(); empty = DStringNew();
	curchromosome = DStringEmtpy();
	result1 = DStringArrayNew(max1+2);
	skip_header(f1,line1,&numfields1,&pos1);
	while (!DStringGetTab(line1,f1,max1,result1,1,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,numfields1,line1,argv[1],&pos1);
		chromosome1 = result1->data+chr1pos;
		sscanf(result1->data[start1pos].string,"%d",&start1);
		sscanf(result1->data[end1pos].string,"%d",&end1);
		NODPRINT("%d\t%s\t%d\t%d",1,Loc_ChrString(chromosome1),start1,end1)
	 	comp = DStringLocCompare(chromosome1,curchromosome);
		if (comp < 0 || (comp == 0 && (start1 < prevstart1 || (start1 == prevstart1 && end1 < prevend1)))) {
			fprintf(stderr,"Cannot annotate because the variant file (%s) is not correctly sorted (sort correctly using \"cg select -s -\")",argv[1]);
			fprintf(stderr,"%s:%d-%d came before %s:%d-%d\n",curchromosome->string,prevstart1,prevend1, chromosome1->string,start1,end1);
			exit(1);
		} else if (comp > 0) {
			DStringCopy(curchromosome,chromosome1);
			nextpos = 0;
		}
		prevstart1 = start1; prevend1 = end1;
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",Loc_ChrString(chromosome1),start1);
			fflush(stderr);
			nextpos += 50000000;
		}
		comp = DStringLocCompare(chromosome2, chromosome1);
		while (comp < 0) {
			if (loaded2) {
				bcol_close(fbcol);
				loaded2 = 0;
			}
			/* keep data of previous */
			/* to avoid allocating new memory everytime, reuse linekeep and associated data */
			chromosomekeep = chromosome2;
			if (argc <= 0) {
				fprintf(stderr,"no bcol file given for chromosome %s\n",chromosome1->string);
				exit(1);
			}
			bcolpos += 2; argc -= 2;
			if (argc < 0) {
				fprintf(stderr,"wrong # of arguments: for each bcol, we must get bcolchr,bcoltype,bcolfile\n");
				exit(1);
			}
			chromosome2 = DStringNew();
			DStringSet(chromosome2,argv[bcolpos]);
			bcolfile = argv[bcolpos+1];
			comp = DStringLocCompare(chromosome2, chromosomekeep);
			DStringDestroy(chromosomekeep);
			if (comp < 0) {
				fprintf(stderr,"Cannot annotate because the bcol file arguments are not correctly sorted");
				fprintf(stderr,"%s came before %s\n",chromosomekeep->string, chromosome2->string);
				exit(1);
			}
			comp = DStringLocCompare(chromosome2, chromosome1);
		}
		if (comp == 0 && !loaded2) {
			fbcol = bcol_open(bcolfile);
			loaded2 = 1;
		}
		if (comp > 0) {
			fprintf(stderr,"no bcol file given for chromosome %s\n",chromosome1->string);
			exit(1);
		}
		if (end1 <= start1) {end1 = start1+1;}
		c = bcol_getbin(fbcol,start1,end1-1);
		if (c < fbcol->typesize) {
			fprintf(stdout,"0\n");
			/* fprintf(stderr,"Could not read value in bcol at position %ld\n",(long int)start1);
			exit(1); */
		} else {
			bcol_printtext(stdout,fbcol->reverse,fbcol->isunsigned,fbcol->type,fbcol->buffer);
			fprintf(stdout,"\n");
		}
	}
	fclose(f1);
	if (loaded2) {
		bcol_close(fbcol);
		loaded2 = 0;
	}
	if (line1) {DStringDestroy(line1);}
	if (linekeep) {DStringDestroy(linekeep);}
	if (empty) {DStringDestroy(empty);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
