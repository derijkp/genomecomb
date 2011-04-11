
#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

int main(int argc, char *argv[]) {
	FILE *f1,*f2 = NULL;
	DString line;
	DString chromosome1,chromosome2,curchromosome;
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	int nchr1=0,start1,end1,nchr2=0,start2,end2;
	int curnchr,curstart,curend;
	int nextpos=0;
	DStringInit(&line);
	DStringInit(&chromosome1);DStringInit(&chromosome2);DStringInit(&curchromosome);
	if ((argc != 9)) {
		fprintf(stderr,"Format is: reg_join file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	if (argv[5][0] != '\0') {
		f2 = fopen64(argv[5],"r");
	}
	chr2pos = atoi(argv[6]);
	start2pos = atoi(argv[7]);
	end2pos = atoi(argv[8]);
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	fprintf(stdout,"chromosome\tbegin\tend\n");
	skip_header(f1, &line);
	get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1);
	if (f2 != NULL) {
		skip_header(f2, &line);
		get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&nchr2,&start2,&end2);
	} else {
		nchr2 = FINISHED;
	}
	if ((nchr1 == FINISHED) && (nchr2 == FINISHED)) exit(0);
	if ((nchr1 < nchr2) || ((nchr1 == nchr2) && (start1 < start2))) {
		DStringCopy(&curchromosome,&chromosome1); curnchr = nchr1; curstart = start1; curend = end1;
		get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1);
	} else {
		DStringCopy(&curchromosome,&chromosome2); curnchr = nchr2; curstart = start2; curend = end2;
		get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&nchr2,&start2,&end2);
	}
	while (1) {
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,chromosome1.string,start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,chromosome2.string,start2,end2);
*/
		if ((nchr1 == FINISHED) && (nchr2 == FINISHED)) break;
		if (nchr1 > curnchr) {
			nextpos = 0;
		}
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",chromosome1.string,start1);
			fflush(stderr);
			nextpos += 10000000;
		}
/*
fprintf(stderr,"# --------------\n");
fprintf(stderr,"# %d\t%d\t%d\n", curnchr,curstart,curend);
fprintf(stderr,"# %d\t%d\t%d\n", nchr1,start1,end1);
fprintf(stderr,"# %d\t%d\t%d\n", nchr2,start2,end2);
*/
		if ((nchr1 < nchr2) || ((nchr1 == nchr2) && (start1 < start2))) {
			if ((nchr1 == curnchr) && (start1 <= curend)) {
				if (end1 > curend) {curend = end1;}
			} else {
				fprintf(stdout,"%s\t%d\t%d\n", curchromosome.string,curstart,curend);
				DStringCopy(&curchromosome,&chromosome1); curnchr = nchr1; curstart = start1; curend = end1;
			}
			get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1);
		} else {
			if ((nchr2 == curnchr) && (start2 <= curend)) {
				if (end2 > curend) {curend = end2;}
			} else {
				fprintf(stdout,"%s\t%d\t%d\n", curchromosome.string,curstart,curend);
				DStringCopy(&curchromosome,&chromosome2); curnchr = nchr2; curstart = start2; curend = end2;
			}
			get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&nchr2,&start2,&end2);
		}
	}
	fprintf(stdout,"%s\t%d\t%d\n", curchromosome.string,curstart,curend);
	fclose(f1);
	if (f2 != NULL) fclose(f2);
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
