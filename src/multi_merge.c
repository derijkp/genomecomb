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
#include "gztools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	VarFile **varfiles = NULL, *varfile;
	Variant *bestvar = NULL;
	DString *alts = NULL;
	int bestpos = 0, bestposnum = 0, *bestposlist = NULL, comp, count, split, i;
	if (argc < 2) {
		fprintf(stderr,"Format is: multi_merge split file1 file2 ...\n");
		exit(EXIT_FAILURE);
	}
	alts = DStringNew();
	split = atoi(argv[1]);
	count = argc - 2;
	varfiles = (VarFile **)malloc(count*sizeof(VarFile *));
	bestposlist = (int *)malloc(count*sizeof(int));
	for (i = 0 ; i < count ; i++) {
		varfiles[i] = OpenVarfile(argv[i+2],split);
		Varfile_checkbasicfields(varfiles[i]);
		varfile_next(varfiles[i]);
	}
	fprintf(stdout,"chromosome\tbegin\tend\ttype\tref\talt\n");
	while (1) {
		/* find best var (first one) */
		bestpos = -1;
		bestposnum = 0;
		for (i = 0 ; i < count ; i++) {
			varfile = varfiles[i];
			if (varfile->error) continue;
			if (bestpos == -1) {
				comp = -1;
			} else {
				comp = varcompare(varfile->var,bestvar,split);
			}
			if (comp == 0) {
				bestposlist[bestposnum++] = i;
			} else if (comp < 0) {
				bestpos = i;
				bestposnum = 0;
				bestposlist[bestposnum++] = i;
				bestvar = varfiles[i]->var;
			}
		}
		if (bestpos == -1) {
			break; /* all files are finished */
		}
		/* print next/best var */
		varputs_chr(bestvar->chr,stdout);
		fprintf(stdout,"\t%d\t%d\t",bestvar->start,bestvar->end);
		DStringputs(bestvar->type,stdout);
		putc_unlocked('\t',stdout);
		DStringputs(bestvar->ref,stdout);
		putc_unlocked('\t',stdout);
		if (split) {
			DStringputs(bestvar->alt,stdout);
		} else {
			DStringCopy(alts,bestvar->alt);
		}
		/* move to next var for all files where current (best) var was found */
		for (i = 0 ; i < bestposnum ; i++) {
			varfile = varfiles[bestposlist[i]];
			bestvar = varfile->var;
			if (!split) {
				/* do we need to add to alts */
				char *alt1 = alts->string,  *alt2 = bestvar->alt->string;
				char *alt1keep,*alt2keep;
				int found;
				alt1keep = alt1;
				while (1) {
					int alt2len;
					alt2keep = alt2;
					while (*alt2 != ',' && *alt2 != '\0') {
						alt2++;
					}
					alt2len = alt2-alt2keep;
					if (alt2len != 1 || *alt2keep != '.') {
						/* check for presence in alt1 */
						alt1keep = alts->string;
						alt1 = alts->string;
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
							if (alts->size == 1 && alts->string[0] == '.') {
								DStringSetS(alts,alt2keep,alt2len);
							} else {
								DStringAppendS(alts,",",1);
								DStringAppendS(alts,alt2keep,alt2len);
							}
						}
					}
					if (*alt2 == '\0') break;
					alt2++;
				}
			}
			varfile_next(varfile);
		}
		if (!split) {
			DStringputs(alts,stdout);
		}
		putc_unlocked('\n',stdout);
	}
	for (i = 0 ; i < count ; i++) {
		CloseVarfile(varfiles[i]);
	}
	exit(EXIT_SUCCESS);
}
