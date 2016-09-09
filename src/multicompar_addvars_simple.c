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
	FILE *orif, *allf, *sregf;
	DStringArray *oriresult=NULL, *allresult=NULL, *sregresult=NULL;
	VariantPos orivarpos, allvarpos, sregvarpos;
	Variant orivar, prevvar, allvar, sregvar;
	DString *oriline = NULL, *allline = NULL, *sregline = NULL;
	DString *ds_q = DStringNew(), *ds_v = DStringNew(), *ds_r = DStringNew(), *ds_u = DStringNew(), *ds_d = DStringNew(), *ds_o = DStringNew();
	DString *out_seq, *out_zyg, *out_a1, *out_a2, *preva1=DStringNew(), *preva2=DStringNew();
	unsigned int *orikeepposs;
	char *orivarsfile, *allvarsfile, *sregfile;
	unsigned int orinumfields,oripos,orierror, numfields;
	unsigned int allnumfields,allpos,allerror;
	unsigned int sregnumfields,sregpos,sregerror;
	int oriseqpos, orizygpos, oria1pos, oria2pos, orikeepsize, orimax = 0, sregmax = 0,prevcomp=-2;
	int split = 1;
	int comp,inregion;
	register int i,j;

	if ((argc < 5)) {
		fprintf(stderr,"Format is: pmulticompar_addvars_simple allvarsfile varsfile split sregfile keepfieldspos1 ...\n");
		exit(EXIT_FAILURE);
	}
	DStringSetS(ds_r,"r",1); DStringSetS(ds_v,"v",1); DStringSetS(ds_u,"u",1);
	DStringSetS(ds_q,"?",1); DStringSetS(ds_d,"-",1); DStringSetS(ds_o,"o",1);
	i = 1;
	allvarsfile = argv[i++];
	orivarsfile = argv[i++];
	split = atoi(argv[i++]);
	sregfile = argv[i++];
	orikeepsize = argc - 5;
	orikeepposs = (unsigned int *)malloc(orikeepsize*sizeof(unsigned int));
	orimax = 0;
	j = 0;
	while (i < argc) {
		orikeepposs[j] = atoi(argv[i]);
		if (orikeepposs[j] > orimax) {orimax = orikeepposs[j];}
		i++; j++;
	}
	/* initialise prevvar */
	prevvar.chr = DStringNew();
	prevvar.type = DStringNew();
	prevvar.ref = DStringNew();
	prevvar.alt = DStringNew();
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
	oria1pos = orivarpos.a1;
	oria2pos = orivarpos.a2;
	if (oriseqpos > orimax) {orimax = oriseqpos;}
	if (orizygpos > orimax) {orimax = orizygpos;}
	if (oria1pos > orimax) {orimax = oria1pos;}
	if (oria2pos > orimax) {orimax = oria2pos;}
	orierror = DStringGetTab(oriline,orif,orimax,oriresult,1,&numfields);
	if (!orierror) {
		check_numfieldserror(numfields,orinumfields,oriline,orivarsfile,&oripos);
		result2var(oriresult,orivarpos,&orivar);
	}
	/* sregfile */
	if (*sregfile != '\0') {
		sregf = fopen64_or_die(sregfile,"r");
		sregline = DStringNew();
		skip_header(sregf,sregline,&sregnumfields,&sregpos);
		sregresult = DStringArrayNew(sregnumfields+2);
		DStringSplitTab(sregline, sregnumfields, sregresult, 0,NULL);
		varpos_fromheader(&sregvarpos,sregresult);
		varpos_max(&(sregvarpos));
		if (sregvarpos.max > sregmax) {sregmax = sregvarpos.max;}
		sregerror = DStringGetTab(sregline,sregf,sregmax,sregresult,1,&numfields);
		if (!sregerror) {
			check_numfieldserror(numfields,sregnumfields,sregline,sregfile,&sregpos);
			result2var(sregresult,sregvarpos,&sregvar);
		}
	} else {
		sregerror = 1;
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
			comp = -2;
		} else {
			comp = varcompare(&allvar,&orivar,split);
		}
		/* comp is -1/1 for same loc/type, but different allele, -2/-2 for different loc/type */
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
			if (oria1pos != -1) {
				putc_unlocked('\t',stdout);
				DStringputs(oriresult->data+oria1pos,stdout);
			}
			if (oria2pos != -1) {
				putc_unlocked('\t',stdout);
				DStringputs(oriresult->data+oria2pos,stdout);
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
			if (split) {
				DStringCopy(prevvar.chr,orivar.chr);
				DStringCopy(prevvar.type,orivar.type);
				/* we do not care about alt differences, so no DStringCopy(prevvar.alt,orivar.alt); */
				prevvar.start = orivar.start;
				prevvar.end = orivar.end;
				DStringCopy(preva1,oriresult->data+oria1pos);
				DStringCopy(preva2,oriresult->data+oria2pos);
			}
			orierror = DStringGetTab(oriline,orif,orimax,oriresult,1,&numfields);
			if (!orierror) {
				check_numfieldserror(numfields,orinumfields,oriline,orivarsfile,&(oripos));
				result2var(oriresult,orivarpos,&orivar);
			}
		} else {
			if (split && comp == -1) {
				/* same location/type, but different allele */
				out_seq = ds_r;
				out_a1 = oriresult->data+oria1pos; out_a2 = oriresult->data+oria2pos;
				if (orizygpos != -1) {
					if (DStringCompare(out_a1,allresult->data+allvarpos.ref) != 0) {
						out_zyg = ds_o;
					} else if (DStringCompare(out_a2,allresult->data+allvarpos.ref) != 0) {
						out_zyg = ds_o;
					} else {
						out_zyg = ds_r;
					}
				}
			} else if (split && prevcomp == 0 && ((comp = varcompare(&allvar,&prevvar,split)) == -1 || comp == 1)) {
				/* previous is from same location, take a1 and a2 from previous */
				out_seq = ds_r;
				out_a1 = preva1; out_a2 = preva2;
				/* keep out_a1 and out_a2 from previous */
				if (orizygpos != -1) {
					if (DStringCompare(out_a1,allresult->data+allvarpos.ref) != 0) {
						out_zyg = ds_o;
					} else if (DStringCompare(out_a2,allresult->data+allvarpos.ref) != 0) {
						out_zyg = ds_o;
					} else {
						out_zyg = ds_r;
					}
				}
				comp = 0;
			} else if (sregline == NULL) {
				out_seq = ds_q; out_zyg = ds_q;
				out_a1 = ds_q; out_a2 = ds_q;
			} else {
				/* search if in sreg sequenced region */
				inregion = 0;
				while (!sregerror) {
					comp = DStringLocCompare(sregvar.chr,allvar.chr);
					if (comp > 0) break;
					if (comp == 0) {
						if (sregvar.end > allvar.start) {
							if (allvar.start >= sregvar.start && allvar.end <= sregvar.end) {
								inregion = 1;
							}
							break;
						}
					}
					sregerror = DStringGetTab(sregline,sregf,sregmax,sregresult,1,&numfields);
					if (!sregerror) {
						check_numfieldserror(numfields,sregnumfields,sregline,sregfile,&(sregpos));
						result2var(sregresult,sregvarpos,&(sregvar));
					}
				}
				if (inregion) {
					out_seq = ds_r; out_zyg = ds_r;
					out_a1 = allresult->data+allvarpos.ref; out_a2 = allresult->data+allvarpos.ref;
				} else {
					out_seq = ds_u; out_zyg = ds_u;
					out_a1 = ds_d; out_a2 = ds_d;
				}
			}
			DStringputs(out_seq,stdout);
			if (orizygpos != -1) {
				putc_unlocked('\t',stdout);
				DStringputs(out_zyg,stdout);
			}
			if (oria1pos != -1) {
				putc_unlocked('\t',stdout);
				DStringputs(out_a1,stdout);
			}
			if (oria2pos != -1) {
				putc_unlocked('\t',stdout);
				DStringputs(out_a2,stdout);
			}
			j = orikeepsize;
			while (j--) {
				putc_unlocked('\t',stdout);
				putc_unlocked('?',stdout);
			}
		}
		prevcomp = comp;
		putc_unlocked('\n',stdout);
		allerror = DStringGetTab(allline,allf,allvarpos.max,allresult,1,&numfields);
		if (allerror) break;
		check_numfieldserror(numfields,allnumfields,allline,allvarsfile,&(allpos));
		result2var(allresult,allvarpos,&allvar);
	}
	exit(EXIT_SUCCESS);
}
