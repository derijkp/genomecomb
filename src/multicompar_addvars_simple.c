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
#include "tools_var.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	FILE *orif, *allf;
	DStringArray *oriresult=NULL, *allresult=NULL;
	VariantPos orivarpos, allvarpos;
	Variant orivar, allvar;
	DString *oriline = NULL, *allline = NULL;
	unsigned int *orikeepposs;
	char *orivarsfile, *allvarsfile;
	unsigned int orinumfields,oripos,orierror, numfields;
	unsigned int allnumfields,allpos,allerror;
	int oriseqpos, orizygpos, orikeepsize, orimax;
	int split = 1;
	int comp;
	register int i,j;

	if ((argc < 4)) {
		fprintf(stderr,"Format is: pmulticompar_addvars_simple allvarsfile varsfile split keepfieldspos1 ...\n");
		exit(EXIT_FAILURE);
	}
	i = 1;
	allvarsfile = argv[i++];
	orivarsfile = argv[i++];
	split = atoi(argv[i++]);
	orikeepsize = argc - 4;
	orikeepposs = (unsigned int *)malloc(orikeepsize*sizeof(unsigned int));
	orimax = 0;
	j = 0;
	while (i < argc) {
		orikeepposs[j] = atoi(argv[i]);
		if (orikeepposs[j] > orimax) {orimax = orikeepposs[j];}
		i++; j++;
	}
	/* orivarfile */
	orif = fopen64_or_die(orivarsfile,"r");
	oriline = DStringNew();
	skip_header(orif,oriline,&orinumfields,&oripos);
	oriresult = DStringArrayNew(orinumfields+2);
	DStringSplitTab(oriline, orinumfields, oriresult, 0,NULL);
	varpos_fromheader(&orivarpos,oriresult);
	varpos_max(&(orivarpos));
	if (orivarpos.max > orimax) {orimax = orivarpos.max;}
	oriseqpos = orivarpos.seq;
	orizygpos = orivarpos.zyg;
	if (oriseqpos > orimax) {orimax = oriseqpos;}
	if (orizygpos > orimax) {orimax = orizygpos;}
	orierror = DStringGetTab(oriline,orif,orimax,oriresult,1,&numfields);
	if (!orierror) {
		check_numfieldserror(numfields,orinumfields,oriline,orivarsfile,&oripos);
		result2var(oriresult,orivarpos,&orivar);
	}
	/* file with all variants */
	varpos_init(&allvarpos);
	allvarpos.chr = 0;
	allvarpos.start = 1;
	allvarpos.end = 2;
	allvarpos.type = 3;
	allvarpos.ref = 4;
	allvarpos.alt = 5;
	varpos_max(&(allvarpos));
	allf = fopen64_or_die(allvarsfile,"r");
	allline = DStringNew();
	allresult = DStringArrayNew(8);
	skip_header(allf,allline,&allnumfields,&allpos);
	allerror = DStringGetTab(allline,allf,6,allresult,1,&numfields);
	if (!allerror) {
		check_numfieldserror(numfields,allnumfields,allline,allvarsfile,&allpos);
		result2var(allresult,allvarpos,&allvar);
	}
	while (!allerror) {
		if (orierror) {
			comp = -1;
		} else {
			comp = varcompare(&allvar,&orivar,split);
		}
		if (comp > 0) {
			fprintf(stderr,"All variants should be in the variant file, so this should not happen");
			exit(EXIT_FAILURE);
		} else if (comp == 0) {
			register unsigned int *cur = orikeepposs;
			if (oriseqpos == -1) {
				putc_unlocked('v',stdout);
			} else {
				DStringputs(oriresult->data+oriseqpos,stdout);
			}
			if (orizygpos != -1) {
				putc_unlocked('\t',stdout);
				DStringputs(oriresult->data+orizygpos,stdout);
			}
			j = orikeepsize;
			while (j--) {
				putc_unlocked('\t',stdout);
				if (*cur != -1) {
					DStringputs(oriresult->data+*cur,stdout);
				} else {
					putc_unlocked('1',stdout);
				}
				cur++;
			}
			orierror = DStringGetTab(oriline,orif,orimax,oriresult,1,&numfields);
			if (!orierror) {
				check_numfieldserror(numfields,orinumfields,oriline,orivarsfile,&(oripos));
				result2var(oriresult,orivarpos,&(orivar));
			}
		} else {
			putc_unlocked('?',stdout);
			if (orizygpos != -1) {
				putc_unlocked('\t',stdout);
				putc_unlocked('?',stdout);
			}
			j = orikeepsize;
			while (j--) {
				putc_unlocked('\t',stdout);
				putc_unlocked('?',stdout);
			}
		}
		putc_unlocked('\n',stdout);
		allerror = DStringGetTab(allline,allf,allvarpos.max,allresult,1,&numfields);
		if (allerror) break;
		check_numfieldserror(numfields,allnumfields,allline,allvarsfile,&(allpos));
		result2var(allresult,allvarpos,&allvar);
	}
	exit(EXIT_SUCCESS);
}
