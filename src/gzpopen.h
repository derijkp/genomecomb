#ifndef GZPOPEN_H_LOADED
#define GZPOPEN_H_LOADED 1

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tools.h"
#include "debug.h"

#define UNCOMPRESSED 1
#define LZ4 2
#define GZ 3
#define RZ 4
#define IN 5
#define ZSTD 6

typedef struct PFILE {
	FILE *f;
	char type;
	DString *filename;
	DString *buffer;
} PFILE;

PFILE *gz_popen(char *filename,char *format);
int gz_pclose(PFILE *file);

#endif
