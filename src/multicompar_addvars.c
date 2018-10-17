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
#include "tools.h"
#include "tools_var.h"
#include "tools_bcol.h"
#include "debug.h"
#include "gztools.h"

typedef struct FieldAnnot {
	char type;
	BCol *bcol;
	VarFile *regf;
	int regfieldpos;
} FieldAnnot;

typedef struct DelInfo {
	DString *chr;
	int start;
	int end;
	int zyg;
} DelInfo;

int inregion(Variant *allvar,VarFile *sreg) {
	int comp;
	while (!sreg->error) {
		comp = DStringLocCompare(sreg->var->chr,allvar->chr);
		if (comp > 0) break;
		if (comp == 0) {
			if (sreg->var->end > allvar->start) {
				if (allvar->start >= sreg->var->start && allvar->end <= sreg->var->end) {
					return 1;
				}
				break;
			}
		}
		varfile_next(sreg,1);
	}
	return 0;
}

void orivarmatch(VarFile *orivars, int *orikeepposs, int orikeepsize, DelInfo *orid) {
	register int *cur = orikeepposs, j;
	int oriseqpos = orivars->varpos.seq;
	int orizygpos = orivars->varpos.zyg;
	int oria1pos = orivars->varpos.a1;
	int oria2pos = orivars->varpos.a2;
	if (oriseqpos == -1) {
		putc_unlocked('v',stdout);
	} else {
		DStringputs(orivars->result->data+oriseqpos,stdout);
	}
	if (orizygpos != -1) {
		putc_unlocked('\t',stdout);
		DStringputs(orivars->result->data+orizygpos,stdout);
	}
	if (oria1pos != -1) {
		putc_unlocked('\t',stdout);
		DStringputs(orivars->result->data+oria1pos,stdout);
	}
	if (oria2pos != -1) {
		putc_unlocked('\t',stdout);
		DStringputs(orivars->result->data+oria2pos,stdout);
	}
	j = orikeepsize;
	while (j--) {
		putc_unlocked('\t',stdout);
		if (*cur != -1) {
			DStringputs(orivars->result->data+*cur,stdout);
		} else {
			putc_unlocked('1',stdout);
		}
		cur++;
	}
	varfile_next(orivars,1);
	if (!orivars->error) {
		if (DStringCompare(orivars->var->chr,orid->chr) != 0) {
			DStringCopy(orid->chr,orivars->var->chr);
			orid->start = -1;
			orid->end = -1;
		}
		if (orivars->var->type->size == 3 && (strncmp(orivars->var->type->string,"del",3) == 0 || strncmp(orivars->var->type->string,"sub",3) == 0)) {
			DStringCopy(orid->chr,orivars->var->chr);
			orid->start = orivars->var->start; orid->end = orivars->var->end;
			if (orizygpos != -1) {
				char c=orivars->result->data[orizygpos].string[0];
				if (c == 'm') {
					orid->zyg = 2;
				} else if (c == 't' || c == 'c') {
					orid->zyg = 1;
				} else {
					orid->zyg = 0;
				}
			}
		}
	}
}

int main(int argc, char *argv[]) {
	VarFile *allvars = NULL;
	VarFile *sreg = NULL;
	VarFile *orivars = NULL;
	VarFile *varallvars = NULL;
	DString *ds_q = DStringNew(), *ds_v = DStringNew(), *ds_r = DStringNew(), *ds_u = DStringNew(), *ds_o = DStringNew(), *ds_c = DStringNew(), *ds_at = DStringNew();
	DString *out_seq=NULL, *out_zyg=NULL, *out_a1=NULL, *out_a2=NULL;
	DelInfo *orid = (DelInfo *)malloc(sizeof(DelInfo));
	FieldAnnot *fieldannotlist = NULL; int checkfieldannot = 0, numbcolannot, numregfiles;
	int *orikeepposs, *varallkeepposs = NULL;
	char *orivarsfile, *allvarsfile, *sregfile, *varallfile;
	int orizygpos, oria1pos, oria2pos, orikeepsize, prevcomp=-2;
	int split = 1;
	int comp, match, pos;
	register int i,j;
	orid->chr = DStringNew();
	orid->start = -1, orid->end = -1, orid->zyg = 2;

	if ((argc < 7)) {
		fprintf(stderr,"Format is: pmulticompar_addvars split allvarsfile varsfile sregfile varallfile numbcolannot ?bcoannotpos bcolaannotfile? ... keepfieldspos1 ...\n");
		exit(EXIT_FAILURE);
	}
	DStringSetS(ds_r,"r",1); DStringSetS(ds_v,"v",1); DStringSetS(ds_u,"u",1);
	DStringSetS(ds_q,"?",1); DStringSetS(ds_o,"o",1);
	DStringSetS(ds_c,"c",1); DStringSetS(ds_at,"@",1);
	i = 1;
	split = atoi(argv[i++]);
	allvarsfile = argv[i++];
	orivarsfile = argv[i++];
	sregfile = argv[i++];
	varallfile = argv[i++];
	numbcolannot = atoi(argv[i++]);
	numregfiles = atoi(argv[i++]);
	orikeepsize = argc - i - 2*numbcolannot - 3*numregfiles;
	fieldannotlist = (FieldAnnot *)malloc(orikeepsize*sizeof(FieldAnnot));
	for (j = 0; j < orikeepsize; j++) {
		fieldannotlist[j].type = '\0';
	}
	for (j = 0; j < numbcolannot; j++) {
		checkfieldannot = 1;
		pos = atoi(argv[i++]);
		fieldannotlist[pos].type = 'b';
		fieldannotlist[pos].bcol = bcol_open(argv[i++]);
	}
	for (j = 0; j < numregfiles; j++) {
		VarFile *regf;
		checkfieldannot = 1;
		pos = atoi(argv[i++]);
		fieldannotlist[pos].type = 'r';
		regf = OpenVarfile(argv[i++],1);
		varfile_next(regf,1);
		fieldannotlist[pos].regf = regf;
		fieldannotlist[pos].regfieldpos = atoi(argv[i++]);
		if (fieldannotlist[pos].regfieldpos > regf->max) {
			regf->max = fieldannotlist[pos].regfieldpos;
		}
	}
	orikeepposs = (int *)malloc(orikeepsize*sizeof(int));
	/* orivarfile */
	orivars = OpenVarfile(orivarsfile,split);
	j = 0;
	while (i < argc) {
		orikeepposs[j] = atoi(argv[i]);
		if (orikeepposs[j] > orivars->max) {orivars->max = orikeepposs[j];}
		i++; j++;
	}
	varfile_next(orivars,1);
	orizygpos = orivars->varpos.zyg;
	oria1pos = orivars->varpos.a1;
	oria2pos = orivars->varpos.a2;
	/* sregfile */
	if (*sregfile != '\0') {
		sreg = OpenVarfile(sregfile,split);
		varfile_next(sreg,1);
	}
	/* varallfile */
	if (*varallfile != '\0') {
		varallvars = OpenVarfile(varallfile,split);
		varallkeepposs = (int *)malloc(orikeepsize*sizeof(int));
		j = 0;
		for (i = 0 ; i < orikeepsize ; i++) {
			DString *el = orivars->header->data+orikeepposs[i];
			varallkeepposs[i] = DStringArraySearch(varallvars->header,el->string,el->size);
			if (varallkeepposs[i] > varallvars->max) {varallvars->max = varallkeepposs[i];}
		}
		/* get first var from varall */
		varfile_next(varallvars,1);
	}
	/* get first variant in ori */
	if (!orivars->error) {
		/* keep end of previous deletion/substitution to check for overlap */
		if (orivars->var->type->size == 3 && (strncmp(orivars->var->type->string,"del",3) == 0 || strncmp(orivars->var->type->string,"sub",3) == 0)) {
			DStringCopy(orid->chr,orivars->var->chr);
			orid->start = orivars->var->start; orid->end = orivars->var->end;
			if (orizygpos != -1) {
				char c=orivars->result->data[orizygpos].string[0];
				if (c == 'm') {
					orid->zyg = 2;
				} else if (c == 't' || c == 'c') {
					orid->zyg = 1;
				} else {
					orid->zyg = 0;
				}
			}
		}
	}
	/* open file with all variants (allvar) */
	allvars = OpenVarfile(allvarsfile,split);
	varfile_next(allvars,1);
	/* go over all vars in allvars, add vars from ori if found, otherwise test sreg, etc. */
	while (!allvars->error) {
		out_seq = ds_q; out_zyg = ds_q; out_a1 = ds_q; out_a2 = ds_q;
		if (orivars->error) {
			comp = -2;
		} else {
			comp = varcompare(allvars->var,orivars->var,split);
		}
		/* comp is -1/1 for same loc/type, but different allele, -2/2 for different loc/type */
		if (comp > 0) {
			fprintf(stderr,"internal error: all variants should be in the (allvar) variant file, so this should not happen");
			exit(EXIT_FAILURE);
		} else if (comp == 0) {
			/* allvar variant found in orivars->var, output this */
			orivarmatch(orivars, orikeepposs, orikeepsize, orid);
		} else if (varallfile[0] != '\0') {
			int refmatch = 0;
			/* var is not in varfile, check in varall for adding data, sreg for u */
			/* find next varall, first pos match for snp, full loc for other */
			while (!varallvars->error) {
				/* find varallvar with same start */
				comp = DStringLocCompare(allvars->var->chr,varallvars->var->chr);
				if (comp == 0) {
					/* if varall hit is a ref, overlap of begin with the ref region is enough */
					/* ref should never overlap a variant */
					if (
						allvars->var->start >= varallvars->var->start
						&& allvars->var->start < varallvars->var->end
						&& varallvars->var->type->size == 3 && varallvars->var->type->string[0] == 'r' && varallvars->var->type->string[1] == 'e' && varallvars->var->type->string[2] == 'f'
					) {
						refmatch = 1;
						break;
					}
					comp = (allvars->var->start - varallvars->var->start);
				}
				if (comp < 0) break;
				if (comp == 0) {
					/* check if we can use the match (chromosome and start) */
					if (varallvars->var->type->size == 3 && varallvars->var->type->string[0] == 's' && varallvars->var->type->string[1] == 'n' && varallvars->var->type->string[2] == 'p') {
						/* for a snp in varall matching only start is enough */
						/* because we use the info at start for indels iun allvars */
						break;
					} else if (allvars->var->end == varallvars->var->end && DStringCompare(allvars->var->type,varallvars->var->type) == 0) {
						/* a full match is also good to stop */
						break;
					}
				}
				varfile_next(varallvars,1);
				if (varallvars->error) {
					comp = -2;
					break;
				}
			}
			if (sreg != NULL) {
				/* search if in sreg sequenced region or not */
				if (inregion(allvars->var,sreg)) {
					out_seq = ds_r;
				} else {
					out_seq = ds_u;
				}
			} else {
				out_seq = ds_q;
			}
			if (refmatch) {out_seq = ds_r;}
			if (comp < 0) {
				/* not in varall, -> unsequenced */
				/* if out_seq was not set based on sreg, set to u */
				match = 0;
				if (out_seq == ds_q) {out_seq = ds_u;}
				if (allvars->var->start < orid->end && allvars->var->end > orid->start && (DStringLocCompare(allvars->var->chr,orid->chr) == 0)) {
					/* overlapping a deletion */
					out_zyg = ds_o;
					if (orid->zyg == 2) {out_a1 = ds_at;} else {out_a1 = allvars->var->ref;}
					if (orid->zyg >= 1) {out_a2 = ds_at;} else {out_a2 = allvars->var->ref;}
				} else if (out_seq == ds_r) {
					out_zyg = ds_r;
					out_a1 = allvars->var->ref;
					out_a2 = allvars->var->ref;
				} else {
					out_zyg = ds_u;
					out_a1 = ds_q; out_a2 = ds_q;
				}
			} else {
				/* comp == 0: hit found in varall, check if it matches completely */
				if (allvars->var->end != varallvars->var->end) {
					match = 1;
				} else if (DStringCompare(allvars->var->type,varallvars->var->type) != 0) {
					match = 1;
				} else if (split && DStringCompare(allvars->var->alt,varallvars->var->alt) != 0) {
					match = 2; /* match except alt allele */
				} else {
					match = 3; /* full match */
				}
				if (match >= 2) {
					/* match of varall, also type, start by using varall seq and zyg (as if match = 3) */
					/* change out_seq later based on data if match == 2 */
					if (varallvars->varpos.seq != -1) {
						/* the value of out_seq i later compared directly to ds_*, so need to set here */
						if (varallvars->result->data[varallvars->varpos.seq].size == 1) {
							if (varallvars->result->data[varallvars->varpos.seq].string[0] == 'u') {
								/* varall is u, so out_seq too */
								out_seq = ds_u;
							} else if (varallvars->result->data[varallvars->varpos.seq].string[0] == 'r') {
								out_seq = ds_r;
							} else if (varallvars->result->data[varallvars->varpos.seq].string[0] == 'v') {
								out_seq = ds_v;
							} else {
								out_seq = varallvars->result->data + varallvars->varpos.seq;
							}
						} else {
							out_seq = varallvars->result->data + varallvars->varpos.seq;
						}
					} else if (out_seq == ds_q) {
						/* if out_seq was not set based on sreg, set to r if varall (different allele) does not have u state */
						out_seq = ds_r;
					}
					if (varallvars->varpos.a1 != -1) {
						out_a1 = varallvars->result->data + varallvars->varpos.a1;
					}
					if (varallvars->varpos.a2 != -1) {
						out_a2 = varallvars->result->data + varallvars->varpos.a2;
					}
					if (match == 3) {
						if (varallvars->varpos.zyg != -1) {
							out_zyg = varallvars->result->data + varallvars->varpos.zyg;
						}
					} else {
						/* partial match in varall, same location/type, but different allele */
						/* match == 2, only for split */
						if (out_seq == ds_u) {
							/* varall is u, so keep out_seq */
						} else if (DStringCompare(out_a1,allvars->result->data + allvars->varpos.alt) != 0
							&& DStringCompare(out_a1,ds_q) != 0
							&& DStringCompare(out_a2,allvars->result->data + allvars->varpos.alt) != 0
							&& DStringCompare(out_a2,ds_q) != 0
							) {
								out_seq = ds_r;
						} else if (out_seq == ds_q) {
							/* if out_seq was not set based on sreg, set to r if varall (different allele) does not have u state */
							out_seq = ds_r;
						}
						if (DStringCompare(out_a1,allvars->result->data+allvars->varpos.ref) != 0) {
							out_zyg = ds_o;
						} else if (DStringCompare(out_a2,allvars->result->data+allvars->varpos.ref) != 0) {
							out_zyg = ds_o;
						} else if (out_zyg == ds_q) {
							out_zyg = ds_r;
						}
						if (out_seq == ds_u)	{out_zyg = ds_u;} /* temp */
					}
				} else {
					/* no match in varall, different end or type  */
					if (out_seq == ds_r) {
						out_zyg = ds_r;
						out_a1 = allvars->var->ref;
						out_a2 = allvars->var->ref;
					} else if (out_seq == ds_u) {
						out_zyg = ds_u; out_a1 = ds_q; out_a2 = ds_q;
					} else {
						out_zyg = ds_u; out_a1 = ds_q; out_a2 = ds_q;
					}
				}
				NODPRINT("%d %d vs %d %d\n",allvars->var->start,allvars->var->end,orid->start,orid->end)
				NODPRINT("ref: %s seq: %s zyg: %s a1: %s a2: %s\n",varallvars->var->ref->string, out_seq->string,out_zyg->string,out_a1->string,out_a2->string)
				if (out_seq == ds_r && allvars->var->start < orid->end && allvars->var->end > orid->start && (DStringLocCompare(allvars->var->chr,orid->chr) == 0)) {
					/* overlapping a deletion */
					if (out_a1 == ds_q || out_a1 == ds_u || DStringCompare(out_a1,varallvars->var->ref) == 0) {
						out_zyg = ds_o;
						if (orid->zyg == 2) {out_a1 = ds_at;} else {out_a1 = allvars->var->ref;}
					}
					if (out_a2 == ds_q || out_a2 == ds_u ||DStringCompare(out_a2,varallvars->var->ref) == 0) {
						out_zyg = ds_o;
						if (orid->zyg == 2 || (out_a1 != ds_at && orid->zyg == 1)) {out_a2 = ds_at;} else {out_a2 = allvars->var->ref;}
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
			if (match) {
				register int *cur = varallkeepposs;
				while (j--) {
					putc_unlocked('\t',stdout);
					if (*cur != -1) {
						DStringputs(varallvars->result->data+*cur,stdout);
					} else {
						putc_unlocked('?',stdout);
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
				out_a1 = orivars->result->data+oria1pos; out_a2 = orivars->result->data+oria2pos;
				if (orizygpos != -1) {
					if (DStringCompare(out_a1,allvars->result->data+allvars->varpos.ref) != 0) {
						out_zyg = ds_o;
					} else if (DStringCompare(out_a2,allvars->result->data+allvars->varpos.ref) != 0) {
						out_zyg = ds_o;
					} else {
						out_zyg = ds_r;
					}
				}
			} else if (split && prevcomp == 0 && ((comp = varcompare(allvars->var,orivars->prevvar,split)) == -1 || comp == 1)) {
				/* previous is from same location/type, take a1 and a2 from previous */
				out_seq = ds_r;
				out_a1 = orivars->prevvar->a1; out_a2 = orivars->prevvar->a2;
				/* keep out_a1 and out_a2 from previous */
				if (orizygpos != -1) {
					if (DStringCompare(out_a1,allvars->result->data+allvars->varpos.ref) != 0) {
						out_zyg = ds_o;
					} else if (DStringCompare(out_a2,allvars->result->data+allvars->varpos.ref) != 0) {
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
				if (inregion(allvars->var,sreg)) {
					out_seq = ds_r; out_zyg = ds_r;
					out_a1 = allvars->result->data+allvars->varpos.ref; out_a2 = allvars->result->data+allvars->varpos.ref;
				} else {
					out_seq = ds_u; out_zyg = ds_u;
					out_a1 = ds_q; out_a2 = ds_q;
				}
			}
			if (out_seq == ds_r && allvars->var->start < orid->end && allvars->var->end > orid->start && (DStringLocCompare(allvars->var->chr,orid->chr) == 0)) {
				/* overlapping a deletion */
				if (out_a1 == ds_q || out_a1 == ds_u || DStringCompare(out_a1,allvars->var->ref) == 0) {
					out_zyg = ds_o;
					if (orid->zyg == 2) {out_a1 = ds_at;} else {out_a1 = allvars->var->ref;}
				}
				if (out_a2 == ds_q || out_a2 == ds_u ||DStringCompare(out_a2,allvars->var->ref) == 0) {
					out_zyg = ds_o;
					if (orid->zyg == 2 || (out_a1 != ds_at && orid->zyg == 1)) {out_a2 = ds_at;} else {out_a2 = allvars->var->ref;}
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
			for (j = 0; j < orikeepsize; j++) {
				if (!checkfieldannot || fieldannotlist[j].type == '\0') {
					putc_unlocked('\t',stdout);
					putc_unlocked('?',stdout);
				} else if (fieldannotlist[j].type == 'b') {
					BCol *bcol = fieldannotlist[j].bcol;
					putc_unlocked('\t',stdout);
					if (bcol_getbinloc(bcol,allvars->var->chr,allvars->var->start,allvars->var->start) != 0) {
						bcol_printtext(stdout,bcol->reverse,bcol->isunsigned,bcol->type,bcol->buffer,bcol->precision);
					} else {
						putc_unlocked('?',stdout);
					}
				} else if (fieldannotlist[j].type == 'r') {
					VarFile *regf = fieldannotlist[j].regf;
					putc_unlocked('\t',stdout);
					if (inregion(allvars->var,regf)) {
						if (fieldannotlist[j].regfieldpos == -1) {
							putc_unlocked('1',stdout);
						} else {
							DStringputs(regf->result->data + fieldannotlist[j].regfieldpos, stdout);
						}
					}
				} else {
					fprintf(stderr, "internal error: wrong fieldannottype\n");
					exit(EXIT_FAILURE);
				}
			}

		}
		prevcomp = comp;
		putc_unlocked('\n',stdout);
		varfile_next(allvars,1);
		if (allvars->error) break;
	}
	CloseVarfile(orivars);
	CloseVarfile(allvars);
	if (varallvars != NULL) CloseVarfile(varallvars);
	exit(EXIT_SUCCESS);
}
