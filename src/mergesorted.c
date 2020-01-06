/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tools.h"
#include "gzpopen.h"
#include "debug.h"

typedef struct MergeFile {
	FILE *f;
	DString *line;
	DStringArray *result;
	int error;
} MergeFile;

int main(int argc, char *argv[]) {
	MergeFile *mergefile, *mergefiles;
	char *header;
	int sortpos[6],numpos,headerline,max,printheader;
	int bestpos = 0, comp, count, i, read;
	char commentchar, *temp;
	if (argc < 2) {
		fprintf(stderr,"Format is: mergesorted commentchar headerline header sortpositions file1 file2 ...\n");
		exit(EXIT_FAILURE);
	}
	commentchar = argv[1][0];
	headerline = atoi(argv[2]);
	header = argv[3];
	temp = argv[4];
	numpos = 0;
	max = 0;
	while (*temp != '\0' && numpos < 6) {
		sortpos[numpos] = atoi(temp);
		if (sortpos[numpos] > max) {max = sortpos[numpos];}
		while (*temp >= 48 && *temp <= 57 && *temp != '\0') temp++;
		while ((*temp <= 48 || *temp >= 57) && *temp != '\0') temp++;
		numpos++;
	}
	if (header[0] != '\0') {
		fprintf(stdout,"%s\n",header);
		printheader = 0;
	} else {
		printheader = 1;
	}
	count = argc - 5;
	argv = argv+5;
	mergefiles = (MergeFile *)malloc(count*sizeof(MergeFile));
	for (i = 0 ; i < count ; i++) {
		FILE *f = gz_popen(argv[i]);
		DString *line=DStringNew();
		mergefile = mergefiles + i;
		mergefile->f = f;
		mergefile->line = line;
		mergefile->result = DStringArrayNew(max+2);
		read = DStringGetLine(line, f);
		while (read != -1) {
			if (line->size > 0 && line->string[0] != commentchar) break;
			if (printheader) {DStringputs(line,stdout); putc_unlocked('\n',stdout);}
			read = DStringGetLine(line, f);
		}
		if (headerline) {
			if (printheader) {DStringputs(line,stdout); putc_unlocked('\n',stdout);}
			read = DStringGetLine(line, f);
		}
		if (printheader) {printheader = 0;}
		if (read == -1) {mergefile->error = 1;} else {mergefile->error = 0;}
		if (!mergefile->error) DStringSplitTab(line,	max, mergefile->result, 0, NULL);
	}
	while (1) {
		MergeFile *bestfile,*testfile;
		DString *beststring, *teststring;
		int filepos;
		/* find best var (first one) */
		bestpos = -1;
		for (filepos = 0 ; filepos < count ; filepos++) {
			mergefile = mergefiles + filepos;
			if (mergefile->error) continue;
			if (bestpos == -1) {
				bestpos = filepos;
				bestfile = mergefiles + bestpos;
			} else {
				testfile = mergefiles + filepos;
				comp = 0;
				for (i = 0 ; i < numpos ; i++) {
					int curpos=sortpos[i];
					beststring = bestfile->result->data+curpos;
					teststring = testfile->result->data+curpos;
					comp = naturalcompare(teststring->string,beststring->string,teststring->size,beststring->size);
					if (comp != 0) break;
				}
				if (comp < 0) {
					bestpos = filepos;
					bestfile = mergefiles + bestpos;
				}
			}
		}
		if (bestpos == -1) {
			break; /* all files are finished */
		}
		/* print next/best var */
		mergefile = mergefiles + bestpos;
		DStringputs(mergefile->line,stdout);
		putc_unlocked('\n',stdout);
		mergefile->error = DStringGetTab(mergefile->line,mergefile->f,max,mergefile->result,0,NULL);
	}
	for (i = 0 ; i < count ; i++) {
		gz_pclose(mergefiles[i].f);
		DStringDestroy(mergefiles[i].line);
		DStringArrayDestroy(mergefiles[i].result);
	}
	exit(EXIT_SUCCESS);
}
