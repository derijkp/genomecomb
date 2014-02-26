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

typedef struct VariantPos {
	int chr;
	int start;
	int end;
	int type;
	int ref;
	int alt;
	int max;
} VariantPos;

typedef struct Variant {
	DString *chr;
	int start;
	int end;
	DString *type;
	DString *ref;
	DString *alt;
} Variant;

/* returns 0 if equal, 1 if equal not including alt, 2 if different */
int varchecksort(Variant *prev,Variant *var,char *filename,int *nextpos) {
	return checksort(prev->chr,&(prev->start),&(prev->end),prev->type,prev->alt,var->chr,var->start,var->end,var->type,var->alt,filename,nextpos);
}

void varputs(Variant var,FILE *f) {
	if (var.chr->size > 3) {
		char *a = var.chr->string;
		if ((a[0] == 'C' || a[0] == 'c') && (a[1] == 'H' || a[1] == 'h') && (a[2] == 'R' || a[2] == 'r')) {
			charputs(a+3,var.chr->size-3,f);
		} else {
			DStringputs(var.chr,f);
		}
	} else {
		DStringputs(var.chr,f);
	}
	fprintf(f,"\t%d\t%d\t",var.start,var.end);
	DStringputs(var.type,f);
	putc_unlocked('\t',f);
	DStringputs(var.ref,f);
	putc_unlocked('\t',f);
	DStringputs(var.alt,f);
	putc_unlocked('\n',f);
}

void result2var(DStringArray *result,VariantPos varpos, Variant *var) {
	var->chr = result->data+varpos.chr;
	sscanf(result->data[varpos.start].string,"%d",&(var->start));
	sscanf(result->data[varpos.end].string,"%d",&(var->end));
	var->type = result->data+varpos.type;
	var->ref = result->data+varpos.ref;
	var->alt = result->data+varpos.alt;
}

void mergealts(
	DString *alt1ds,DString *alt2ds
) {
	char *alt1 = alt1ds->string,  *alt2 = alt2ds->string;
	char *alt1keep,*alt2keep;
	int found;
	alt1keep = alt1;
	while (1) {
		alt2keep = alt2;
		while (*alt2 != ',' && *alt2 != '\0') {
			alt2++;
		}
		/* check for presence in alt1 */
		alt1keep = alt1ds->string;
		alt1 = alt1ds->string;
		found = 0;
		while (1) { 
			while (*alt1 != ',' && *alt1 != '\0') {
				alt1++;
			}
			if ((alt1-alt1keep == alt2-alt2keep) && (strncmp(alt1keep,alt2keep,alt1-alt1keep) == 0)) {
				found = 1;
				break;
			}
			if (*alt1 == '\0') break;
			alt1++;
			alt1keep = alt1;
		}
		if (!found) {
			DStringAppendS(alt1ds,",",1);
			DStringAppendS(alt1ds,alt2keep,alt2-alt2keep);
		}
		if (*alt2 == '\0') break;
		alt2++;
	}
}

int main(int argc, char *argv[]) {
	FILE *f1,*f2;
	DStringArray *result1=NULL,*result2=NULL;
	VariantPos var1pos,var2pos;
	Variant prev1,prev2;
	Variant var1,var2;
	DString *line1 = NULL,*line2 = NULL;
	char *filename1,*filename2;
	unsigned int numfields1,numfields2,numfields,pos1,pos2;
	int split = 1;
	int comp,comptype,compalt = 0,compcheck;
	int error2,nextpos=0,same = 0;
	prev1.chr = DStringNew(); prev1.start = -1; prev1.end = -1; prev1.type = DStringNew(); prev1.ref = DStringNew(); prev1.alt = DStringNew();
	prev2.chr = DStringNew(); prev2.start = -1; prev2.end = -1; prev2.type = DStringNew(); prev2.ref = DStringNew(); prev2.alt = DStringNew();
	var1.ref = NULL; var2.ref = NULL;
	var1.chr=NULL,var2.chr=NULL,var1.type = NULL,var2.type = NULL,var1.alt = NULL,var2.alt = NULL;

	if ((argc != 16)) {
		fprintf(stderr,"Format is: multi_merge file1 chrpos1 startpos1 endpos1 type1pos ref1pos alt1pos file2 chrpos2 startpos2 endpos2 type2pos ref2pos alt2pos split");
		exit(EXIT_FAILURE);
	}
	filename1 = argv[1];
	f1 = fopen64_or_die(filename1,"r");
	var1pos.chr = atoi(argv[2]);
	var1pos.start = atoi(argv[3]);
	var1pos.end = atoi(argv[4]);
	var1pos.type = atoi(argv[5]);
	var1pos.ref = atoi(argv[6]);
	var1pos.alt = atoi(argv[7]);
	var1pos.max = var1pos.chr;
	if (var1pos.start > var1pos.max) {var1pos.max = var1pos.start;} ; if (var1pos.end > var1pos.max) {var1pos.max = var1pos.end;} ;
	if (var1pos.type > var1pos.max) {var1pos.max = var1pos.type;} ;
	if (var1pos.ref > var1pos.max) {var1pos.max = var1pos.ref;} ; 
	if (var1pos.alt > var1pos.max) {var1pos.max = var1pos.alt;} ;
	line1 = DStringNew(); line2=DStringNew();
	/* The following allocation is not destroyed at end as it may point to something else */
	/* This will leak mem, but as the prog is finished anyway ... */
	result1 = DStringArrayNew(var1pos.max+2);
	filename2 = argv[8];
	f2 = fopen64_or_die(filename2,"r");
	var2pos.chr = atoi(argv[9]);
	var2pos.start = atoi(argv[10]);
	var2pos.end = atoi(argv[11]);
	var2pos.type = atoi(argv[12]);
	var2pos.ref = atoi(argv[13]);
	var2pos.alt = atoi(argv[14]);
	split = atoi(argv[15]);
NODPRINT("var_annot %s %d %d %d %d %d %s %d %d %d %d %d ...",
	filename1,var1pos.chr,var1pos.start,var1pos.end,var1pos.type,var1pos.alt,
	filename2,var2pos.chr,var2pos.start,var2pos.end,var2pos.type,var2pos.alt
);
	var2pos.max = var2pos.chr ; if (var2pos.start > var2pos.max) {var2pos.max = var2pos.start;} ; if (var2pos.end > var2pos.max) {var2pos.max = var2pos.end;} ;
	if (var2pos.type > var2pos.max) {var2pos.max = var2pos.type;} ;
	if (var2pos.ref > var2pos.max) {var2pos.max = var2pos.ref;} ; 
	if (var2pos.alt > var2pos.max) {var2pos.max = var2pos.alt;} ;
	result2 = DStringArrayNew(var2pos.max+2);
	skip_header(f1,line1,&numfields1,&pos1);
	skip_header(f2,line2,&numfields2,&pos2);
	fprintf(stdout,"chromosome\tbegin\tend\ttype\tref\talt\n");
	error2 = DStringGetTab(line2,f2,var2pos.max,result2,1,&numfields); pos2++;
	if (!error2) {
		check_numfieldserror(numfields,numfields2,line2,filename2,&pos2);
		result2var(result2,var2pos,&var2);
		varchecksort(&prev2,&var2,filename2,&nextpos);
	}
NODPRINT("line2 %s,%d,%d %s",Loc_ChrString(chromosome2),start2,end2,line2->string)
	while (!DStringGetTab(line1,f1,var1pos.max,result1,1,&numfields)) {
		pos1++;
		check_numfieldserror(numfields,numfields1,line1,filename1,&pos1);
		result2var(result1,var1pos,&var1);
		NODPRINT("line1 (a=%3.3s) %s,%d,%d %s",type1->string,chromosome1->string,start1,end1,alt1->string)
		/*
		fprintf(stdout,"----- %d\t%s\t%d\t%d\n",1,Loc_ChrString(chromosome1),start1,end1);
		fprintf(stdout,"--------- %d\t%s\t%d\t%d\n",2,Loc_ChrString(chromosome2),start2,end2);
		*/
		compcheck = varchecksort(&prev1,&var1,filename1,&nextpos);
		if (!split && compcheck < 2) {
			fprintf(stderr,"error in \"%s\": file uses split alleles (\"%s %d %d %s\" occurs more than once and you are not running multicompar with the -split option)",
				filename1,prev1.chr->string,prev1.start,prev1.end,prev1.type->string);
			exit(EXIT_FAILURE);
		}
		/* if (var1.start >= nextpos) {
			fprintf(stderr, "%s-%d\n",Loc_ChrString(var1.chr),var1.start);
			fflush(stderr);
			nextpos += 50000000;
		} */
		while (!error2) {
			same = 0;
			NODPRINT("line2 %s,%d,%d %s %s",Loc_ChrString(prevchromosome2),prevstart2,prevend2,prevtype2->string,prevalt2->string)
			comp = DStringLocCompare(var2.chr,var1.chr);
			if (comp == 0) {
				if (var2.start == var1.start) {
					if (var2.end == var1.end) {
						comptype = DStringCompare(var2.type, var1.type);
						if (comptype == 0) {
							if (split) {
								compalt = DStringCompare(var2.alt, var1.alt);
								if (compalt == 0) {
									same = 1; break;
								} else if (compalt > 0) break;
							} else {
								same = 1; break;
							}
						} else if (comptype > 0) break; 
					} else if (var2.end > var1.end) break; 
				} else if (var2.start > var1.start) break;
			} else if (comp > 0) break;
			varputs(var2,stdout);
			error2 = DStringGetTab(line2,f2,var2pos.max,result2,1,&numfields); pos2++;
			if (error2)  {
				same = 0;
				break;
			} else {
				check_numfieldserror(numfields,numfields2,line2,filename2,&pos2);
			}
			result2var(result2,var2pos,&var2);
			compcheck = varchecksort(&prev2,&var2,filename2,&nextpos);
			if (!split && compcheck < 2) {
				fprintf(stderr,"error in \"%s\": file uses split alleles (\"%s %d %d %s\" occurs more than once and you are not running multicompar with the -split option)",
					filename2,prev2.chr->string,prev2.start,prev2.end,prev2.type->string);
				exit(EXIT_FAILURE);
			}
		}
		if (!error2 && same) {
			/* merge alts */
			if (!split) {
				mergealts(var1.alt,var2.alt);
			}
			error2 = DStringGetTab(line2,f2,var2pos.max,result2,1,&numfields); pos2++;
			if (!error2) {
				check_numfieldserror(numfields,numfields2,line2,filename2,&pos2);
				result2var(result2,var2pos,&var2);
				compcheck = varchecksort(&prev2,&var2,filename2,&nextpos);
				if (!split && compcheck < 2) {
					fprintf(stderr,"error in \"%s\": file uses split alleles (\"%s %d %d %s\" occurs more than once and you are not running multicompar with the -split option)",
						filename2,prev2.chr->string,prev2.start,prev2.end,prev2.type->string);
					exit(EXIT_FAILURE);
				}
			}
		}
		varputs(var1,stdout);
	}
	while (!error2) {
		varputs(var2,stdout);
		error2 = DStringGetTab(line2,f2,var2pos.max,result2,1,&numfields); pos2++;
		if (error2)  {
			break;
		} else {
			check_numfieldserror(numfields,numfields2,line2,filename2,&pos2);
		}
		result2var(result2,var2pos,&var2);
		compcheck = varchecksort(&prev2,&var2,filename2,&nextpos);
		if (!split && compcheck < 2) {
			fprintf(stderr,"error in \"%s\": file uses split alleles (\"%s %d %d %s\" occurs more than once and you are not running multicompar with the -split option)",
				filename2,prev2.chr->string,prev2.start,prev2.end,prev2.type->string);
			exit(EXIT_FAILURE);
		}
	}
	fclose(f1);
	fclose(f2);
	DStringDestroy(prev1.chr); DStringDestroy(prev1.type); DStringDestroy(prev1.alt);
	DStringDestroy(prev2.chr); DStringDestroy(prev2.type); DStringDestroy(prev2.alt);
	if (line1) {DStringDestroy(line1);}
	if (line2) {DStringDestroy(line2);}
	if (result1) {DStringArrayDestroy(result1);}
	if (result2) {DStringArrayDestroy(result2);}
	exit(EXIT_SUCCESS);
}
