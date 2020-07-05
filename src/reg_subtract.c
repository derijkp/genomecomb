/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "tools.h"

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	DString line;
	DString *chromosome1 = NULL,*chromosome2 = NULL,curchromosome;
	DString *prevchromosome1 = DStringNew(),*prevchromosome2 = DStringNew();
	int chr1pos,start1pos,end1pos,chr2pos,start2pos,end2pos,max1,max2,comp;
	int start1,end1,start2,end2;
	int prevstart1 = -1, prevend1 = -1, prevstart2 = -1,prevend2 = -1;
	int error,error2;
#ifdef SHOWPROGRESS
	int nextpos=0;
#endif
	DStringInit(&line);
	chromosome1 = DStringNew();
	chromosome2 = DStringNew();
	DStringInit(&curchromosome);
	if ((argc != 9)) {
		fprintf(stderr,"Format is: reg_subtract file1 chrpos1 startpos1 endpos1 file2 chrpos2 startpos2 endpos2");
		exit(EXIT_FAILURE);
	}
	f1 = fopen64_or_die(argv[1],"r");
	chr1pos = atoi(argv[2]);
	start1pos = atoi(argv[3]);
	end1pos = atoi(argv[4]);
	max1 = chr1pos ; if (start1pos > max1) {max1 = start1pos;} ; if (end1pos > max1) {max1 = end1pos;} ;
	f2 = fopen64_or_die(argv[5],"r");
	chr2pos = atoi(argv[6]);
	start2pos = atoi(argv[7]);
	end2pos = atoi(argv[8]);
NODPRINT("reg_subtract %s %d %d %d %s %d %d %d",argv[1],chr1pos,start1pos,end1pos,argv[5],chr2pos,start2pos,end2pos)
	max2 = chr2pos ; if (start2pos > max2) {max2 = start2pos;} ; if (end2pos > max2) {max2 = end2pos;} ;
	fprintf(stdout,"chromosome\tbegin\tend\n");
	skip_header(f2,&line,NULL,NULL);
	skip_header(f1,&line,NULL,NULL);
	error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1);
	if (error) {exit(0);}
	error2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&start2,&end2);
	while (1) {
/*
fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,chromosome1,start1,end1);
fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,chromosome2,start2,end2);
*/
	 	comp = DStringLocCompare(chromosome1,&curchromosome);
		if (comp > 0) {
			DStringCopy(&curchromosome,chromosome1);
#ifdef SHOWPROGRESS
			nextpos = 0;
#endif
		}
#ifdef SHOWPROGRESS
		if (chromosome1 != NULL && start1 >= nextpos) {
			fprintf(stderr, "%s-%d\n",Loc_ChrString(chromosome1),start1);
			fflush(stderr);
			nextpos += 10000000;
		}
#endif
	 	comp = DStringLocCompare(chromosome2,chromosome1);
		if ((comp < 0) || ((comp == 0) && (end2 <= start1))) {
			error2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&start2,&end2);
			if (error2)  {
				fprintf(stdout,"%s\t%d\t%d\n", Loc_ChrString(chromosome1),start1,end1);
				break;
			}
			checksortreg(prevchromosome2,&prevstart2,&prevend2,chromosome2,start2,start2,argv[5]);
		} else if ((comp > 0) || ((comp == 0) && (end1 <= start2))) {
			fprintf(stdout,"%s\t%d\t%d\n", Loc_ChrString(chromosome1),start1,end1);
			error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1);
			if (error) break;
			checksortreg(prevchromosome1,&prevstart1,&prevend1,chromosome1,start1,start1,argv[1]);
		} else {
			if (start2 > start1) {
				fprintf(stdout,"%s\t%d\t%d\n", Loc_ChrString(chromosome1),start1,start2);
			}
			if (end2 >= end1) {
				error = get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1);
				if (error) break;
				checksortreg(prevchromosome1,&prevstart1,&prevend1,chromosome1,start1,start1,argv[1]);
			} else {
				start1 = end2;
				error2 = get_region(f2,&line,chr2pos,start2pos,end2pos,max2,&chromosome2,&start2,&end2);
				if (error2)  {
					fprintf(stdout,"%s\t%d\t%d\n", Loc_ChrString(chromosome1),start1,end1);
					break;
				}
				checksortreg(prevchromosome2,&prevstart2,&prevend2,chromosome2,start2,start2,argv[5]);
			}
		}
	}
	while (!get_region(f1,&line,chr1pos,start1pos,end1pos,max1,&chromosome1,&start1,&end1)) {
		checksortreg(prevchromosome1,&prevstart1,&prevend1,chromosome1,start1,start1,argv[1]);
		fprintf(stdout,"%s\t%d\t%d\n", Loc_ChrString(chromosome1),start1,end1);
	}
	FCLOSE(f1);
	FCLOSE(f2);
	DStringClear(&line);
	exit(EXIT_SUCCESS);
}
