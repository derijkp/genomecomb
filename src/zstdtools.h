#ifndef ZSTDTOOLS_H_LOADED
#define ZSTDTOOLS_H_LOADED 1

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#define ZSTD_STATIC_LINKING_ONLY
#include "zstd.h"
#include "debug.h"

typedef  uint8_t BYTE;

#define KB *(1<<10)
#define MB *(1<<20)
#define GB *(1<<30)

#define _1BIT  0x01
#define _2BITS 0x03
#define _3BITS 0x07
#define _4BITS 0x0F
#define _8BITS 0xFF

#define MAGICNUMBER_SIZE    4
#define ZSTDIO_MAGICNUMBER   0x184D2204
#define ZSTDIO_SKIPPABLE0    0x184D2A50
#define ZSTDIO_SKIPPABLEMASK 0xFFFFFFF0
#define BUFFER_EMPTY UINT_MAX
#define ZSTDv06_MAGICNUMBER	0xFD2FB526   /* v0.6 */
#define ZSTDv07_MAGICNUMBER	0xFD2FB527   /* v0.7 */

typedef struct zstdres {
	FILE *finput;
	unsigned version;
	size_t compressedsize;
	unsigned long long contentsize;
	unsigned long long framepos;
	unsigned long long int inframepos;
	unsigned long long framefilepos;
	unsigned int frameread;
	char *buffer;
	unsigned int buffersize;
	unsigned long long int bufferpos;
	int lastblock;
	char *outbuffer;
	unsigned int outbuffersize;
	unsigned int currentblock;
	uint64_t currentpos;
	unsigned long long framesize;
	FILE *findex;
	uint64_t indexstart;
	uint64_t indexbsize;
	uint64_t usize;
	ZSTD_frameHeader zfh;
} ZSTDres;

#define ERROR(message) {\
	fprintf(stderr,"%s\n",message);\
	fflush(stderr);\
	exit(1);\
}
#define READ(f,size,var,errormsg) {\
	if (fread(var, 1, size, f) != size) ERROR(errormsg);\
}

char *zstd_findindex(char *filename);
unsigned ZSTDIO_readLE32 (const void* s);
ZSTDres *zstdopen(FILE* finput,FILE *indexfile);
int zstd_readheader(ZSTDres *res);
int zstd_readframe(ZSTDres *res);
int zstd_skipframe(ZSTDres *res);
void zstdclose(ZSTDres *res);
void zstdindex_write(FILE *findex,off_t pos);
uint64_t zstdindex_read(FILE *findex);

ZSTDres *zstd_openfile(char *file);
int zstd_seek(ZSTDres *res, uint64_t pos, int where);
int zstd_read(ZSTDres *res, void *data, uint64_t size);
int zstd_get(ZSTDres *res);
#define zstd_close(res) zstdclose(res)

#endif
