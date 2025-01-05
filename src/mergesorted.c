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
	PFILE *pf;
	DString *line;
	DString *sort;
	DStringArray *result;
	int error;
	struct MergeFile *next;
} MergeFile;

void makesortfield(MergeFile *mergefile,int sortpos[6],int numpos) {
	DString *sort = mergefile->sort;
	DString *field;
	char *a;
	int i;
	int curpos=sortpos[0],alen;
	field = mergefile->result->data+curpos;
	a = field->string;
	alen = field->size;
	if (alen >= 3) {
		if ((a[0] == 'C' || a[0] == 'c') && (a[1] == 'H' || a[1] == 'h') && (a[2] == 'R' || a[2] == 'r')) {
			a += 3; alen -= 3;
			if (alen && a[0] == '-') {
				a++; alen--;
			}
		}
	}
	DStringSetS(sort,a,alen);
	for (i = 1 ; i < numpos ; i++) {
		curpos=sortpos[i];
		field = mergefile->result->data+curpos;
		DStringAppendS(sort," ",1);
		DStringAppendS(sort,field->string,field->size);
	}
}


MergeFile *insertitem(MergeFile *firstmergefile,MergeFile *insertmergefile) {
	MergeFile *testfile, *prev;
	DString *insertstring, *teststring;
	int comp;
	if (firstmergefile == NULL) {
		insertmergefile->next = NULL;
		return insertmergefile;
	}
	insertstring = insertmergefile->sort;
	prev = NULL;
	testfile = firstmergefile;
	while (1) {
		teststring = testfile->sort;
		comp = naturalcompare(insertstring->string,teststring->string,insertstring->size,teststring->size);
		if (comp <= 0) {
			insertmergefile->next = testfile;
			if (prev == NULL) {
				return insertmergefile;
			} else {
				prev->next = insertmergefile;
				return firstmergefile;
			}
		}
		if (testfile->next == NULL) {
			break;
 		}
		prev = testfile;
		testfile = testfile->next;
	}
	insertmergefile->next = NULL;
	testfile->next = insertmergefile;
	return firstmergefile;
}

int main(int argc, char *argv[]) {
	MergeFile *mergefile, *mergefiles, *firstmergefile;
	char *header;
	int sortpos[6],numpos,headerline,max,printheader;
	int count, i, read;
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
	firstmergefile = NULL;
	for (i = 0 ; i < count ; i++) {
		PFILE *pf = gz_popen(argv[i],NULL);
		DString *line = DStringNew();
		mergefile = mergefiles + i;
		mergefile->pf = pf;
		mergefile->line = line;
		mergefile->sort = DStringNew();
		mergefile->result = DStringArrayNew(max+2);
		read = DStringGetLine(line, pf->f);
		while (read != -1) {
			if (line->size > 0 && line->string[0] != commentchar) break;
			if (printheader) {DStringputs(line,stdout); putc_unlocked('\n',stdout);}
			read = DStringGetLine(line, pf->f);
		}
		if (headerline) {
			if (printheader) {DStringputs(line,stdout); putc_unlocked('\n',stdout);}
			read = DStringGetLine(line, pf->f);
		}
		if (printheader) {printheader = 0;}
		if (read == -1) {mergefile->error = 1;} else {mergefile->error = 0;}
		if (!mergefile->error) {
			DStringSplitTab(line,	max, mergefile->result, 0, NULL);
			makesortfield(mergefile,sortpos,numpos);
			firstmergefile = insertitem(firstmergefile,mergefile);
		}
	}
	while (firstmergefile != NULL) {
		/* print top var */
		mergefile = firstmergefile;
		DStringputs(mergefile->line,stdout);
		putc_unlocked('\n',stdout);
		mergefile->error = DStringGetTab(mergefile->line,mergefile->pf->f,max,mergefile->result,0,NULL);
		if (mergefile->error) {
			if (mergefile->next == NULL) break;
			firstmergefile = mergefile->next;
		} else {
			firstmergefile = firstmergefile->next;
			makesortfield(mergefile,sortpos,numpos);
			firstmergefile = insertitem(firstmergefile,mergefile);
		}
	}
	for (i = 0 ; i < count ; i++) {
		gz_pclose(mergefiles[i].pf);
		DStringDestroy(mergefiles[i].line);
		DStringDestroy(mergefiles[i].sort);
		DStringArrayDestroy(mergefiles[i].result);
	}
	exit(EXIT_SUCCESS);
}
