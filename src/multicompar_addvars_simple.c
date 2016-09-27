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
#include "gztools.h"

typedef struct VarFile {
	char *file;
	GZFILE *f;
	DString *line;
	DStringArray *result;
	VariantPos varpos;
	Variant var;
	int max;
	unsigned int numfields;
	unsigned int pos;
	unsigned int error;
} VarFile;

VarFile *OpenVarfile(char *filename) {
	unsigned int numfields;
	VarFile *varfile = malloc(sizeof(VarFile));
	varfile->file = filename;
	varfile->f = gz_open(filename);
	varfile->line = DStringNew();
	gz_skip_header(varfile->f,varfile->line,&varfile->numfields,&varfile->pos);
	varfile->result = DStringArrayNew(varfile->numfields+2);
	DStringSplitTab(varfile->line, varfile->numfields, varfile->result, 0,NULL);
	varpos_fromheader(&varfile->varpos,varfile->result);
	varpos_max(&(varfile->varpos));
	varfile->max = varfile->varpos.max;
	varfile->error = gz_DStringGetTab(varfile->line,varfile->f,varfile->max,varfile->result,1,&numfields);
	if (!varfile->error) {
		check_numfieldserror(numfields,varfile->numfields,varfile->line,varfile->file,&varfile->pos);
		result2var(varfile->result,varfile->varpos,&varfile->var);
	}
	return varfile;
}

int inregion(Variant *allvar,VarFile *sreg) {
	unsigned int numfields;
	int comp;
	while (!sreg->error) {
		comp = DStringLocCompare(sreg->var.chr,allvar->chr);
		if (comp > 0) break;
		if (comp == 0) {
			if (sreg->var.end > allvar->start) {
				if (allvar->start >= sreg->var.start && allvar->end <= sreg->var.end) {
					return 1;
				}
				break;
			}
		}
		sreg->error = gz_DStringGetTab(sreg->line,sreg->f,sreg->max,sreg->result,1,&numfields);
		if (!sreg->error) {
			check_numfieldserror(numfields,sreg->numfields,sreg->line,sreg->file,&(sreg->pos));
			result2var(sreg->result,sreg->varpos,&(sreg->var));
		}
	}
	return 0;
}

int main(int argc, char *argv[]) {
	VarFile *sreg = NULL;
	GZFILE *orif=NULL, *allf=NULL, *varallf=NULL;
	DStringArray *oriresult=NULL, *allresult=NULL, *varallresult=NULL;
	VariantPos orivarpos, allvarpos, varallvarpos;
	Variant orivar, prevvar, allvar, varallvar;
	DString *oriline = NULL, *allline = NULL, *varallline = NULL;
	DString *ds_q = DStringNew(), *ds_v = DStringNew(), *ds_r = DStringNew(), *ds_u = DStringNew(), *ds_d = DStringNew(), *ds_o = DStringNew(), *ds_c = DStringNew(), *ds_at = DStringNew();
	DString *out_seq, *out_zyg, *out_a1, *out_a2, *preva1=DStringNew(), *preva2=DStringNew(), *oridchr = DStringNew();
	DString *allref = NULL;
	unsigned int *orikeepposs, *varallkeepposs;
	char *orivarsfile, *allvarsfile, *sregfile, *varallfile;
	int oridstart = -1, oridend = -1, oridzyg = 2;
	unsigned int orinumfields,oripos,orierror, numfields;
	unsigned int allnumfields,allpos,allerror;
	unsigned int varallnumfields,varallpos,varallerror;
	int oriseqpos, orizygpos, oria1pos, oria2pos, orikeepsize, orimax = 0, varallmax = 0,prevcomp=-2;
	int split = 1;
	int comp;
	register int i,j;

	if ((argc < 6)) {
		fprintf(stderr,"Format is: pmulticompar_addvars_simple allvarsfile varsfile split sregfile varallfile keepfieldspos1 ...\n");
		exit(EXIT_FAILURE);
	}
	DStringSetS(ds_r,"r",1); DStringSetS(ds_v,"v",1); DStringSetS(ds_u,"u",1);
	DStringSetS(ds_q,"?",1); DStringSetS(ds_d,"-",1); DStringSetS(ds_o,"o",1);
	DStringSetS(ds_c,"c",1); DStringSetS(ds_at,"@",1);
	i = 1;
	allvarsfile = argv[i++];
	orivarsfile = argv[i++];
	split = atoi(argv[i++]);
	sregfile = argv[i++];
	varallfile = argv[i++];
	orikeepsize = argc - i;
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
	orif = gz_open(orivarsfile);
	oriline = DStringNew();
	gz_skip_header(orif,oriline,&orinumfields,&oripos);
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
	/* sregfile */
	if (*sregfile != '\0') {
		sreg = OpenVarfile(sregfile);
	}
	/* varallfile */
	if (*varallfile != '\0') {
		varallf = gz_open(varallfile);
		varallline = DStringNew();
		gz_skip_header(varallf,varallline,&varallnumfields,&varallpos);
		varallresult = DStringArrayNew(varallnumfields+2);
		DStringSplitTab(varallline, varallnumfields, varallresult, 0,NULL);
		varpos_fromheader(&varallvarpos,varallresult);
		varpos_max(&(varallvarpos));
		if (varallvarpos.max > varallmax) {varallmax = varallvarpos.max;}
		varallkeepposs = (unsigned int *)malloc(orikeepsize*sizeof(unsigned int));
		varallmax = 0;
		j = 0;
		for (i = 0 ; i < orikeepsize ; i++) {
			DString *el = oriresult->data+orikeepposs[i];
			varallkeepposs[i] = DStringArraySearch(varallresult,el->string,el->size);
			if (orikeepposs[i] > varallmax) {varallmax = orikeepposs[i];}
		}
		/* get first var from varall */
		varallerror = gz_DStringGetTab(varallline,varallf,varallmax,varallresult,1,&numfields);
		if (!varallerror) {
			check_numfieldserror(numfields,varallnumfields,varallline,varallfile,&varallpos);
			result2var(varallresult,varallvarpos,&varallvar);
		}
	} else {
		varallerror = 1;
	}
	/* get first variant in ori */
	orierror = gz_DStringGetTab(oriline,orif,orimax,oriresult,1,&numfields);
	if (!orierror) {
		check_numfieldserror(numfields,orinumfields,oriline,orivarsfile,&oripos);
		result2var(oriresult,orivarpos,&orivar);
		/* keep end of previous deletion/substitution to check for overlap */
		if (orivar.type->size == 3 && (strncmp(orivar.type->string,"del",3) == 0 || strncmp(orivar.type->string,"sub",3) == 0)) {
			DStringCopy(oridchr,orivar.chr);
			oridstart = orivar.start; oridend = orivar.end;
			if (orizygpos != -1) {
				char c=oriresult->data[orizygpos].string[0];
				if (c == 'm') {
					oridzyg = 2;
				} else if (c == 't' || c == 'c') {
					oridzyg = 1;
				} else {
					oridzyg = 0;
				}
			}
		}
	}
	/* open file with all variants (allvar) */
	varpos_init(&allvarpos);
	allvarpos.chr = 0;
	allvarpos.start = 1;
	allvarpos.end = 2;
	allvarpos.type = 3;
	allvarpos.ref = 4;
	allvarpos.alt = 5;
	varpos_max(&(allvarpos));
	allf = gz_open(allvarsfile);
	allline = DStringNew();
	allresult = DStringArrayNew(8);
	gz_skip_header(allf,allline,&allnumfields,&allpos);
	allerror = gz_DStringGetTab(allline,allf,6,allresult,1,&numfields);
	if (!allerror) {
		check_numfieldserror(numfields,allnumfields,allline,allvarsfile,&allpos);
		result2var(allresult,allvarpos,&allvar);
	}
	/* go over all vars in allvar, add vars from ori if found, otherwise test sreg, etc. */
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
			/* allvar variant found in orivar, output this */
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
			orierror = gz_DStringGetTab(oriline,orif,orimax,oriresult,1,&numfields);
			if (!orierror) {
				check_numfieldserror(numfields,orinumfields,oriline,orivarsfile,&(oripos));
				result2var(oriresult,orivarpos,&orivar);
				if (DStringCompare(orivar.chr,oridchr) != 0) {
					DStringCopy(oridchr,orivar.chr);
					oridend = -1;
				}
				if (orivar.type->size == 3 && (strncmp(orivar.type->string,"del",3) == 0 || strncmp(orivar.type->string,"sub",3) == 0)) {
					DStringCopy(oridchr,orivar.chr);
					oridstart = orivar.start; oridend = orivar.end;
					if (orizygpos != -1) {
						char c=oriresult->data[orizygpos].string[0];
						if (c == 'm') {
							oridzyg = 2;
						} else if (c == 't' || c == 'c') {
							oridzyg = 1;
						} else {
							oridzyg = 0;
						}
					}
				}
			}
		} else if (varallline != NULL) {
			/* var is not in varfile, check in varall for adding data, sreg for u */
			while (!varallerror) {
				/* find varallvar with same start */
				comp = DStringLocCompare(allvar.chr,varallvar.chr);
				if (comp == 0) {
					comp = allvar.start - varallvar.start;
				}
				if (varallvar.type->size != 3 || varallvar.type->string[0] != 's') {
					/* for something other than a SNP in the varall, we need a full location match to stop */
					if (comp <= 0 && allvar.end == varallvar.end && DStringCompare(allvar.type,varallvar.type) == 0) break;
				} else {
					/* use only start match if snp */
					if (comp <= 0) break;
				}
				varallerror = gz_DStringGetTab(varallline,varallf,varallmax,varallresult,1,&numfields);
				if (!varallerror) {
					check_numfieldserror(numfields,varallnumfields,varallline,varallfile,&(varallpos));
					result2var(varallresult,varallvarpos,&varallvar);
				} else {
					comp = -2;
				}
			}
			if (sreg != NULL) {
				/* search if in sreg sequenced region or not */
				if (inregion(&allvar,sreg)) {
					out_seq = ds_r;
				} else {
					out_seq = ds_u;
				}
			} else {
				out_seq = ds_q;
			}
			allref = allresult->data + allvarpos.ref;
			if (comp != 0) {
				/* not in varall, -> unsequenced */
				/* if out_seq was not set based on sreg, set to u */
				if (out_seq == ds_q) {out_seq = ds_u;}
				if (allvar.end <= oridend && allvar.start >= oridstart) {
					/* overlapping a deletion */
					out_zyg = ds_o;
					if (oridzyg == 2) {out_a1 = ds_at;} else {out_a1 = allref;}
					if (oridzyg >= 1) {out_a2 = ds_at;} else {out_a2 = allref;}
				} else if (out_seq == ds_r) {
					out_zyg = ds_r;
					out_a1 = allref;
					out_a2 = allref;
				} else {
					out_zyg = ds_u;
					out_a1 = ds_d; out_a2 = ds_d;
				}
			} else {
				if (allvar.end != varallvar.end) {comp = 2;}
				if (comp == 0) {
					if (DStringCompare(allvar.type,varallvar.type) != 0) {comp = 2;}
					if (comp == 0 && split) {
						comp = DStringCompare(allvar.alt,varallvar.alt);
						if (comp != 0) {comp = 1;}
					} else {comp = 1;}
				}
				if (comp == 0) {
					/* complete match of varall, also allele, so use varall seq and zyg */
					out_seq = varallresult->data + varallvarpos.seq;
					out_zyg = varallresult->data + varallvarpos.zyg;
					out_a1 = varallresult->data + varallvarpos.a1;
					out_a2 = varallresult->data + varallvarpos.a2;
				} else if (comp == 1) {
					/* partial match in varall, same location/type, but different allele */
					out_a1 = varallresult->data+varallvarpos.a1;
					out_a2 = varallresult->data+varallvarpos.a2;
					if (varallvarpos.seq != -1 && varallresult->data[varallvarpos.seq].size == 1 && varallresult->data[varallvarpos.seq].string[0] == 'u') {
						/* varall is u, so out_seq too */
						out_seq = ds_u;
					} else if (DStringCompare(out_a1,allresult->data + allvarpos.alt) != 0
						&& DStringCompare(out_a1,ds_q) != 0
						&& DStringCompare(out_a1,ds_d) != 0
						&& DStringCompare(out_a2,allresult->data + allvarpos.alt) != 0
						&& DStringCompare(out_a2,ds_q) != 0
						&& DStringCompare(out_a2,ds_d) != 0
						) {
							out_seq = ds_r;
					} else if (out_seq == ds_q) {
						/* if out_seq was not set based on sreg, set to r if varall (different allele) does not have u state */
						out_seq = ds_r;
					}
					if (orizygpos != -1) {
						if (DStringCompare(out_a1,allresult->data+allvarpos.ref) != 0) {
							out_zyg = ds_o;
						} else if (DStringCompare(out_a2,allresult->data+allvarpos.ref) != 0) {
							out_zyg = ds_o;
						} else {
							out_zyg = ds_r;
						}
					}
				} else {
					/* no match in varall, different end or type  */
					if (out_seq == ds_r) {
						out_zyg = ds_r;
						out_a1 = allref;
						out_a2 = allref;
					} else if (out_seq == ds_u) {
						out_zyg = ds_u; out_a1 = ds_d; out_a2 = ds_d;
					} else {
						out_zyg = ds_u; out_a1 = ds_q; out_a2 = ds_q;
					}
				}
NODPRINT("%d %d vs %d %d\n",allvar.start,allvar.end,oridstart,oridend)
NODPRINT("ref: %s seq: %s zyg: %s a1: %s a2: %s\n",varallvar.ref->string, out_seq->string,out_zyg->string,out_a1->string,out_a2->string)
				if (out_seq == ds_r && allvar.end <= oridend && allvar.start >= oridstart) {
					/* overlapping a deletion */
					if (out_a1 == ds_q || out_a1 == ds_u || DStringCompare(out_a1,varallvar.ref) == 0) {
						out_zyg = ds_o;
						if (oridzyg == 2) {out_a1 = ds_at;} else {out_a1 = allref;}
					}
					if (out_a2 == ds_q || out_a2 == ds_u ||DStringCompare(out_a2,varallvar.ref) == 0) {
						out_zyg = ds_o;
						if (oridzyg == 2 || (out_a1 != ds_at && oridzyg == 1)) {out_a2 = ds_at;} else {out_a2 = allref;}
					}
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
			if (comp > -2) {
				register unsigned int *cur = varallkeepposs;
				while (j--) {
					putc_unlocked('\t',stdout);
					if (*cur != -1) {
						DStringputs(varallresult->data+*cur,stdout);
					} else {
						putc_unlocked('1',stdout);
					}
					cur++;
				}
			} else {
				while (j--) {
					putc_unlocked('\t',stdout);
					putc_unlocked('?',stdout);
				}
			}
		} else {
			/* no match in orivar and no varall present: */
			/* check for same loc/diff allele, and if not found check sreg (if present) */
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
				/* previous is from same location/type, take a1 and a2 from previous */
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
			} else if (sreg == NULL) {
				/* no sreg file available, all to unknown */
				out_seq = ds_q; out_zyg = ds_q;
				out_a1 = ds_q; out_a2 = ds_q;
			} else {
				/* search if in sreg sequenced region or not */
				if (inregion(&allvar,sreg)) {
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
		allerror = gz_DStringGetTab(allline,allf,allvarpos.max,allresult,1,&numfields);
		if (allerror) break;
		check_numfieldserror(numfields,allnumfields,allline,allvarsfile,&(allpos));
		result2var(allresult,allvarpos,&allvar);
	}
	gz_close(orif);
	gz_close(allf);
	gz_close(varallf);
	exit(EXIT_SUCCESS);
}
