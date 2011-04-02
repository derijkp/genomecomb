
#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	DString line;
	DString chromosome1,chromosome2;
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2;
	int nchr1=0,start1,end1,nchr2=0,start2,end2;
	int error,error2,curchr=0,nextpos=0;
	DStringInit(&line);DStringInit(&chromosome1);DStringInit(&chromosome2);
	if ((argc != 9)) {
		fprintf(stderr,"Format is: reg_subtract file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2");
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
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	fprintf(stdout,"chromosome\tbegin\tend\n");
	DStringGetLine(&line,f2);
	DStringGetLine(&line, f1);
	error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1);
	if (error) {exit(0);}
	error2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&nchr2,&start2,&end2);
	while (1) {
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,chromosome1,start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,chromosome2,start2,end2);
*/
		if (nchr1 > curchr) {
			curchr = nchr1;
			nextpos = 0;
		}
		if (start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",chromosome1.string,start1);
			fflush(stderr);
			nextpos += 10000000;
		}
		if ((nchr2 < nchr1) || ((nchr2 == nchr1) && (end2 <= start1))) {
			error2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&nchr2,&start2,&end2);
			if (error2)  {
				fprintf(stdout,"%s\t%d\t%d\n", chromosome1.string,start1,end1);
				break;
			}
		} else if ((nchr1 < nchr2) || ((nchr1 == nchr2) && (end1 <= start2))) {
			fprintf(stdout,"%s\t%d\t%d\n", chromosome1.string,start1,end1);
			error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1);
			if (error) break;
		} else {
			if (start2 > start1) {
				fprintf(stdout,"%s\t%d\t%d\n", chromosome1.string,start1,start2);
			}
			if (end2 >= end1) {
				error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1);
				if (error) break;
			} else {
				start1 = end2;
				error2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&nchr2,&start2,&end2);
				if (error2)  {
					fprintf(stdout,"%s\t%d\t%d\n", chromosome1.string,start1,end1);
					break;
				}
			}
		}
	}
	while (!get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&nchr1,&start1,&end1)) {
		fprintf(stdout,"%s\t%d\t%d\n", chromosome1.string,start1,end1);
	}
	fclose(f1);
	fclose(f2);
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
