#ifndef GZPOPEN_LOADED
#define GZPOPEN_LOADED 1

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <ctype.h>
#include "debug.h"
#include "tools.h"
#include "gzpopen.h"

FILE *gz_popen(char *filename) {
	int type = UNCOMPRESSED, len = strlen(filename);
	DString *buffer=DStringNew();
	FILE *result;
	char *cur;
	if (len == 1 && *filename == '-') {
		type = IN;
	} else if (len >= 3) {
		cur = filename + len - 1;
		if (len >= 4 && *cur == '4') {
			if (*(--cur) == 'z' && *(--cur) == 'l' && *(--cur) == '.') {
				type = LZ4;
			}
		} else if (len >= 4 && *cur == 't') {
			if (*(--cur) == 's' && *(--cur) == 'z' && *(--cur) == '.') {
				type = ZSTD;
			}
		} else if (*cur == 'z') {
			cur--;
			if (*(cur) == 'g' && *(--cur) == '.') {
				type = GZ;
			} else if (*(cur) == 'r' && *(--cur) == '.') {
				type = RZ;
			}
		}
	}
	if (type == IN) {
		return stdin;
	} else if (type == UNCOMPRESSED) {
		return fopen64_or_die(filename,"r");
	} else if (type == LZ4) {
		DStringSet(buffer,"lz4 -q -d -c ");
		DStringAppend(buffer,filename);
		return popen(buffer->string,"r");
	} else if (type == ZSTD) {
		DStringSet(buffer,"zstd-mt -k -q -d -c ");
		DStringAppend(buffer,filename);
		return popen(buffer->string,"r");
	} else if (type == RZ || type == GZ) {
		DStringSet(buffer,"zcat ");
		DStringAppend(buffer,filename);
		return popen(buffer->string,"r");
	} else {
		return fopen64_or_die(filename,"r");
	}
	return result;
}

int gz_pclose(FILE *file) {
	int r;
	r = pclose(file);
	if (r == -1) {
		if (errno == ECHILD) {
			return fclose(file);
		} else {
			fprintf(stderr, "%s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
	return r;
}

#endif
