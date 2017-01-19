#ifndef LZ4TOOLS_H_LOADED
#define LZ4TOOLS_H_LOADED 1

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "lz4.h"
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
#define LZ4IO_MAGICNUMBER   0x184D2204
#define LZ4IO_SKIPPABLE0    0x184D2A50
#define LZ4IO_SKIPPABLEMASK 0xFFFFFFF0
#define LEGACY_MAGICNUMBER  0x184C2102
#define BUFFER_EMPTY UINT_MAX

typedef struct lz4res {
	FILE *finput;
	unsigned version;
	unsigned blockMode;
	unsigned blockChecksumFlag;
	unsigned contentSizeFlag;
	unsigned contentChecksumFlag;
	unsigned int blocksize;
	char *buffer;
	char *outbuffer;
	unsigned int currentblock;
	uint64_t currentpos;
	char *writebuffer;
	unsigned int writesize;
	FILE *findex;
	unsigned int startblocks;
	uint64_t indexstart;
	uint64_t indexbsize;
	uint64_t usize;
} LZ4res;

#define ERROR(message) {\
	fprintf(stderr,"%s\n",message);\
	fflush(stderr);\
	exit(1);\
}
#define READ(f,size,var,errormsg) {\
	if (fread(var, 1, size, f) != size) ERROR(errormsg);\
}

char *lz4_findindex(char *filename);
unsigned LZ4IO_readLE32 (const void* s);
LZ4res *lz4open(FILE* finput,FILE *indexfile);
int lz4_readblock(LZ4res *res, unsigned int startblock);
void lz4close(LZ4res *res);
void lz4index_write(FILE *findex,off_t pos);
uint64_t lz4index_read(FILE *findex);

LZ4res *lz4_openfile(char *file);
void lz4_seek(LZ4res *res, int64_t pos, int where);
int lz4_read(LZ4res *res, void *data, uint64_t size);
int lz4_get(LZ4res *res);
#define lz4_close(res) lz4close(res)

#endif
