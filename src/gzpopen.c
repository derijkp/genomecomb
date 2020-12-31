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

PFILE *gz_popen(char *filename,char *format) {
	int type = UNCOMPRESSED, len = strlen(filename);
	PFILE *result = (PFILE *)malloc(sizeof(PFILE));
	DString *buffer = NULL;
	char *cur;
	result->filename = DStringNewFromChar(filename);
	result->buffer = DStringNew();
	buffer = result->buffer;
	if (len == 1 && *filename == '-') {
		result->f = stdin;
		result->type = 's';
		return result;
	} else if (len >= 3) {
		char *keepcur = filename + len - 1;
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
		if (format != NULL) {
			int flen = strlen(format);
			int elen = 0;
			if (flen == 3 && strncmp(format,"sam",3) == 0) {
				while (*cur != '.' && len > 0) {
					cur--; len--; elen++;
				}
				cur++;
				if (*cur == '.') {len = 0;}
				if (
					(elen == 3 && strncmp(cur,"bam",elen) == 0)
					|| (elen == 4 && strncmp(cur,"cram",elen) == 0)
				) {
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
			result->f = fopen64_or_die(filename,"r");
			result->type = 'f';
		} else {
			result->f = popen(buffer->string,"r");
			result->type = 'p';
		}
	}
	if (result->f == NULL) {
		fprintf(stderr,"Could not open file %s using popen (%s): %s\n",filename,buffer->string,strerror(errno));
		exit(1);
	}
	return result;
}

int gz_pclose(PFILE *pfile) {
	int r;
	if (pfile->type == 'p') {
		r = pclose(pfile->f);
	} else if (pfile->type == 'f') {
		r = fclose(pfile->f);
	}
	if (r == -1) {
		if (errno == ECHILD) {
			FCLOSE(pfile->f);
		} else {
			if (pfile->type == 'p') {
				fprintf(stderr, "error closing stream \"%s\": %s\n", pfile->buffer->string, strerror(errno));
			} else if (pfile->type == 'f') {
				fprintf(stderr, "error closing file %s: %s\n", pfile->filename->string, strerror(errno));
			} else {
			}
			fprintf(stderr, "%s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
	DStringDestroy(pfile->buffer);
	DStringDestroy(pfile->filename);
	free(pfile);
	return r;
}

#endif
