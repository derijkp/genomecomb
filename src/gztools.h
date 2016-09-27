#ifndef GZTOOLS_H_LOADED
#define GZTOOLS_H_LOADED 1

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tools.h"
#include "lz4tools.h"
#include "razf.h"
#include "debug.h"

#define UNCOMPRESSED 1
#define LZ4 2
#define GZ 3
#define RZ 4
#define STDIN 5

typedef struct GZFILE {
	char *filename;
	int type;
	FILE *fun;
	LZ4res *lz4;
	RAZF *rz;
} GZFILE;

GZFILE *gz_open(char *filename);
void gz_seek(GZFILE *f, int64_t pos, int where);
int gz_read(GZFILE *f,void *data,uint64_t size);
int gz_get(GZFILE *f);
void gz_close(GZFILE *f);

int gz_DStringGetLine(DString *dstring,	GZFILE *f1);
void gz_skip_header(GZFILE *f1, DString *linePtr,unsigned int *numfields,unsigned int *pos);
int gz_DStringGetTab(DString *line,	GZFILE *f1, int max, DStringArray *result, int setzero,unsigned int *numfields);

#endif
