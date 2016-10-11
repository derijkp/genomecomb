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
#include "debug.h"
#include "tools.h"

int main(int argc, char *argv[]) {
	FILE *f1,*f2 = NULL;
	DString line;
	DString *chromosome1 = NULL,*chromosome2 = NULL,curchromosome;
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	int start1,end1,start2,end2,finished1 = 0, finished2 = 0,comp;
	int curstart,curend,cursamechr;
#ifdef SHOWPROGRESS
	int nextpos=0;
#endif
	DStringInit(&line);
	chromosome1 = DStringNew();
	chromosome2 = DStringNew();
	DStringInit(&curchromosome);
	if ((argc != 9)) {
		fprintf(stderr,"Format is: reg_join file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64_or_die(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	if (argv[5][0] != '\0') {
		f2 = fopen64_or_die(argv[5],"r");
	}
	chr2pos = atoi(argv[6]);
	start2pos = atoi(argv[7]);
	end2pos = atoi(argv[8]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
NODPRINT("reg_join %s %d %d %d %s %d %d %d",argv[1],chr1pos,start1pos,end1pos,argv[5],chr2pos,start2pos,end2pos)
	fprintf(stdout,"chromosome\tbegin\tend\n");
	skip_header(f1, &line,NULL,NULL);
	finished1 = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1);
	if (f2 != NULL) {
		skip_header(f2, &line,NULL,NULL);
		finished2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&start2,&end2);
	} else {
		finished2 = 1;
		DStringDestroy(chromosome2);
		chromosome2 = NULL;
	}
	if (finished1 && finished2) exit(0);
	comp = DStringLocCompare(chromosome1,chromosome2);
	if ((comp < 0) || ((comp == 0) && (start1 < start2))) {
		DStringCopy(&curchromosome,chromosome1); curstart = start1; curend = end1;
		finished1 = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1);
	} else {
		DStringCopy(&curchromosome,chromosome2); curstart = start2; curend = end2;
		finished2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&start2,&end2);
	}
	while (1) {
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,Loc_ChrString(chromosome1),start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,Loc_ChrString(chromosome2),start2,end2);
*/
		if (finished1 && finished2) break;
	 	comp = DStringLocCompare(chromosome1,&curchromosome);
		if (comp > 0) {
#ifdef SHOWPROGRESS
			int nextpos=0;
#endif
			cursamechr = 0;
		} else if (comp == 0) {
			cursamechr = 1;
		} else {
			cursamechr = 0;
		}
#ifdef SHOWPROGRESS
		if (!finished1 && start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",Loc_ChrString(chromosome1),start1);
			fflush(stderr);
			nextpos += 10000000;
		}
#endif
/*
fprintf(stderr,"# --------------\n");
fprintf(stderr,"# %s\t%d\t%d\n", curchromosome.string,curstart,curend);
fprintf(stderr,"# %s\t%d\t%d\n", Loc_ChrString(chromosome1),start1,end1);
fprintf(stderr,"# %s\t%d\t%d\n", Loc_ChrString(chromosome2),start2,end2);
*/
	 	comp = DStringLocCompare(chromosome1,chromosome2);
		if ((comp < 0) || ((comp == 0) && (start1 < start2))) {
			if ((cursamechr) && (start1 <= curend)) {
				if (end1 > curend) {curend = end1;}
			} else {
				fprintf(stdout,"%s\t%d\t%d\n", curchromosome.string,curstart,curend);
				DStringCopy(&curchromosome,chromosome1); curstart = start1; curend = end1;
			}
			finished1 = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1);
		 	comp = DStringLocCompare(chromosome1,&curchromosome);
			if (comp < 0 || (comp == 0 && start1 < curstart)) {
				fprintf(stderr,"file1 is not correctly sorted (sort correctly using \"cg select -s -\")\n");
				exit(1);
			}
		} else {
		 	comp = DStringLocCompare(chromosome2,&curchromosome);
			if ((comp == 0) && (start2 <= curend)) {
				if (end2 > curend) {curend = end2;}
			} else {
				fprintf(stdout,"%s\t%d\t%d\n", curchromosome.string,curstart,curend);
				DStringCopy(&curchromosome,chromosome2); curstart = start2; curend = end2;
			}
			finished2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&start2,&end2);
		 	comp = DStringLocCompare(chromosome2,&curchromosome);
			if (comp < 0 || (comp == 0 && start2 < curstart)) {
				fprintf(stderr,"file2 is not correctly sorted (sort correctly using \"cg select -s -\")\n");
				exit(1);
			}
		}
	}
	fprintf(stdout,"%s\t%d\t%d\n", curchromosome.string,curstart,curend);
	fclose(f1);
	if (f2 != NULL) fclose(f2);
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
