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
#include "gzpopen.h"

FILE *sam_popen(char *filename) {
	int len = strlen(filename);
	DString *buffer=DStringNew();
	FILE *result;
	char *cur,*keepcur;
	if (len == 1 && *filename == '-') {
		return stdin;
	} else if (len >= 3) {
		keepcur = filename + len - 1;
		cur = keepcur;
		if (len >= 4 && *cur == '4') {
			if (*(--cur) == 'z' && *(--cur) == 'l' && *(--cur) == '.') {
				DStringAppend(buffer,"lz4 -q -d -c ");
				DStringAppend(buffer,filename);
				keepcur = --cur; len -= 4;
			} else {
				cur = keepcur;
			}
		} else if (len >= 4 && *cur == 't') {
			if (*(--cur) == 's' && *(--cur) == 'z' && *(--cur) == '.') {
				DStringAppend(buffer,"zstd-mt -T 1 -k -q -d -c ");
				DStringAppend(buffer,filename);
				keepcur = --cur; len -= 4;
			} else {
				cur = keepcur;
			}
		} else if (*cur == 'z') {
			cur--;
			if ((*(cur) == 'g' || *(cur) == 'r') && *(--cur) == '.') {
				DStringAppend(buffer,"zcat ");
				DStringAppend(buffer,filename);
				keepcur = --cur; len -= 3;
			} else {
				cur = keepcur;
			}
		}
		if (len >= 4 && *cur == 'm') {
			if (*(--cur) == 'a' && *(--cur) == 'b' && *(--cur) == '.') {
				if (buffer->size == 0) {
					DStringAppend(buffer,"samtools view --no-PG -h ");
					DStringAppend(buffer,filename);
				} else {
					DStringAppend(buffer,"| samtools view --no-PG -h ");
				}
				keepcur = --cur;
			} else {
				cur = keepcur;
			}
		}
	}
	if (buffer->size == 0) {
		return fopen64_or_die(filename,"r");
	} else {
		return popen(buffer->string,"r");
	}
	return result;
}

int sam_pclose(FILE *file) {
	int r;
	r = pclose(file);
	if (r == -1) {
		if (errno == ECHILD) {
			FCLOSE(file);
		} else {
			fprintf(stderr, "%s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
	return r;
}

int main(int argc, char *argv[]) {
	FILE *f = NULL;
	int pos;
	register int c, first = 1;
	if ((argc <= 1)) {
		fprintf(stderr,"Format is: samcat ?-header header? samfile1 ...\n");
		exit(EXIT_FAILURE);
	}
	if (strlen(argv[1]) == 7 && strncmp(argv[1],"-header",7) == 0) {
		if (strlen(argv[2]) > 0) {
			fprintf(stdout,"%s\n",argv[2]);
		}
		pos = 3;
	} else {
		f = sam_popen(argv[1]);
		while ((c=getc_unlocked(f))!=EOF) {
			putc_unlocked(c,stdout);
		}
		gz_pclose(f);
		pos = 2;
	}
	while (pos < argc) {
		f = sam_popen(argv[pos++]);
		while ((c=getc_unlocked(f))!=EOF) {
			if (first) {
				if (c != '@') break;
				first = 0;
			}
			if (c == '\n') {first = 1;}
		}
		if (c != EOF) putc_unlocked(c,stdout);
		while ((c=getc_unlocked(f))!=EOF) {
			putc_unlocked(c,stdout);
		}
		sam_pclose(f);
	}
	exit(EXIT_SUCCESS);
}
