/* #define DEBUG 1 */

#include "zstdtools.h"

/* already defined in tools.c */
/*
FILE *fopen64_or_die(char *filename,char *mode) {
	FILE *f;
	f = fopen64(filename,mode);
	if (f == NULL) {
		fprintf(stderr,"Error opening file %s: %s.\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return(f);
}
*/
/* already defined in tools.h, but we do not load tools.h here */
FILE *fopen64_or_die(char *filename,char *mode);

uint32_t ZSTDIO_readLE24 (const void* s) {
	const unsigned char* srcPtr = (const unsigned char*)s;
	unsigned value32 = srcPtr[0];
	value32 += (srcPtr[1]<<8);
	value32 += (srcPtr[2]<<16);
	return value32;
}

uint32_t ZSTDIO_readLE32 (const void* s) {
	const unsigned char* srcPtr = (const unsigned char*)s;
	unsigned value32 = srcPtr[0];
	value32 += (srcPtr[1]<<8);
	value32 += (srcPtr[2]<<16);
	value32 += ((unsigned)srcPtr[3])<<24;
	return value32;
}

char *zstd_findindex(char *filename) {
	char *outfile = (char *)malloc(strlen(filename) + 7);
	sprintf(outfile,"%s.zsti",filename);
	return outfile;
}

#define yaml_checkparam(buffer,read,start,test,size) { \
if (read < (start+size) || strncmp(buffer+start,test,size) != 0 || (buffer[start+size] != ' ' && buffer[start+size] != '\n')) { \
	fprintf(stderr,"zstdi error: only %s supported\n",test); \
	exit(1); \
} \
}

#ifdef DEBUG
void debug_printprogress(ZSTDres *res,char *msg) {
	DPRINT("* %s at %zu  currentpos: %zu framepos: %llu contentsize:%llu inframepos:%llu type:%d",msg,ftello(res->finput),res->currentpos,res->framepos,res->contentsize,res->inframepos,res->zfh.frameType);
}
#endif

#define MIN_ZSTD_HEADERSIZE ZSTD_FRAMEHEADERSIZE_MIN(ZSTD_f_zstd1)
// #define MIN_ZSTD_HEADERSIZE ZSTD_FRAMEHEADERSIZE_MIN

/* returns contentsize, which is 0 for skippedframes */
int zstd_readheader(ZSTDres *res) {
	size_t error;
	DPRINT("zstd_readheader ftel start:%ld",ftell(res->finput));
	if (res->inframepos != 0 && res->frameread == 0) {
		fprintf(stderr, "read or skip frame before reading header of new one\n");
		exit(1);
	}
	res->framepos = res->framepos + res->contentsize;
	res->frameread = 0;
	res->framefilepos = ftell(res->finput);
	error = fread(res->buffer, 1, MIN_ZSTD_HEADERSIZE, res->finput);
	DPRINT("tel2:%ld",ftell(res->finput));
	DPRINT("error:%lu",error);
	DPRINT("MIN_ZSTD_HEADERSIZE:%d",MIN_ZSTD_HEADERSIZE);
	if (error == 0 && feof(res->finput)) return(0);
	if (error != MIN_ZSTD_HEADERSIZE) ERROR("error in file: incomplete zstd frame header");
	error = ZSTD_getFrameHeader(&(res->zfh), (const void*) res->buffer, MIN_ZSTD_HEADERSIZE);
	DPRINT("ZSTD_getFrameHeader:%lu",error);
	if (ZSTD_isError(error)) {
		fprintf(stderr,"%s\n",ZSTD_getErrorName(error));
		exit(1);
	} else if (error > MIN_ZSTD_HEADERSIZE) {
		res->inframepos = error;
		error = fread(res->buffer+MIN_ZSTD_HEADERSIZE, 1, res->inframepos-MIN_ZSTD_HEADERSIZE, res->finput);
		DPRINT("tel3:%ld",ftell(res->finput));
		DPRINT("fread:%ld",error);
		if (error != res->inframepos-MIN_ZSTD_HEADERSIZE) ERROR("could not read zstd header");
		error = ZSTD_getFrameHeader(&(res->zfh), (const void*) res->buffer, res->inframepos);
	} else {
		res->inframepos = MIN_ZSTD_HEADERSIZE;
	}
	DPRINT("tel4:%ld",ftell(res->finput));
	if (ZSTD_isError(error)) {
		fprintf(stderr,"error reading zstd frame header: %s\n",ZSTD_getErrorName(error));
		exit(1);
	}
	res->contentsize = ZSTD_getFrameContentSize((const void *)res->buffer, res->inframepos);
	if (res->contentsize == ZSTD_CONTENTSIZE_UNKNOWN) {
		fprintf(stderr,"error cannot determine size of zstd frame, use zstdmt with -T option to compress\n");
		exit(1);
	} else if (res->contentsize == ZSTD_CONTENTSIZE_ERROR) {
		fprintf(stderr,"error determining size of zstd frame\n");
		exit(1);
	}
#ifdef DEBUG
printf("contentsize:%lld\n",res->contentsize);
printf("telend:%ld\n",ftell(res->finput));
debug_printprogress(res,"readheader");
#endif
	return(res->contentsize);
}

uint32_t zstd_readblockheader(ZSTDres *res) {
	uint32_t cBlockHeader;
	uint32_t blocksize;
	char *blockbuffer;
	if (feof(res->finput)) return(0);
	DPRINT("---- readblockheader at ftello: %zu ----",ftello(res->finput));
	if (res->bufferpos + 3 > res->buffersize) {
		res->buffersize += 3;
		res->buffer = (char *)realloc(res->buffer,res->buffersize);
	}
	blockbuffer = res->buffer + res->bufferpos;
	READ(res->finput,3,blockbuffer,"error in file: incomplete zstd block header");
	cBlockHeader = ZSTDIO_readLE24(blockbuffer);
       blocksize = cBlockHeader >> 3;
       res->lastblock = cBlockHeader & 1;
	return(blocksize);
}

int zstd_readblock(ZSTDres *res, size_t blocksize) {
	if (feof(res->finput)) return(0);
	if (res->bufferpos + blocksize + 3 > res->buffersize) {
		res->buffersize += blocksize + 3;
		res->buffer = (char *)realloc(res->buffer,res->buffersize);
	}
	if (res->lastblock && res->zfh.checksumFlag) {
		blocksize += 4;
	}
	READ(res->finput,blocksize,res->buffer + res->bufferpos + 3,"error in file: incomplete zstd block");
	res->bufferpos += 3 + blocksize;
	return(1);
}

int zstd_readframe(ZSTDres *res) {
	size_t numdecompressed;
	uint32_t blocksize;
#ifdef DEBUG
debug_printprogress(res,"readframe");
#endif
	if (res->frameread) return(1);
	if (res->zfh.frameType == ZSTD_skippableFrame) {
		fprintf(stderr,"error: trying to read skippable zstd frame at %jd",ftello(res->finput));
	}
	res->bufferpos = res->inframepos;
	res->lastblock = 0;
	while (1) {
		blocksize = zstd_readblockheader(res);
		if (!blocksize) break;
		zstd_readblock(res,blocksize);
		if (res->lastblock) break;
	}
	res->inframepos = 0;
	res->frameread = 1;
	if (res->outbuffersize < res->contentsize) {
		res->outbuffer = realloc(res->outbuffer,res->contentsize);
		res->outbuffersize = res->contentsize;
	}
	numdecompressed = ZSTD_decompress(res->outbuffer, res->contentsize, res->buffer, res->bufferpos);
	if (ZSTD_isError(numdecompressed)) {
		fprintf(stderr,"error uncompressing zstd: %s\n",ZSTD_getErrorName(numdecompressed));
		exit(1);
	}
	NODPRINT("numdecompressed: %d",numdecompressed);
	return(1);
}

#if defined (__cplusplus) || (defined (__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) /* C99 */)
# include <stdint.h>
	typedef uint32_t U32;
#else
	typedef unsigned int        U32;
#endif

unsigned MEM_isLittleEndian(void)
{
	const union { U32 u; BYTE c[4]; } one = { 1 };   /* don't use static : performance detrimental  */
	return one.c[0];
}

U32 MEM_read32(const void* memPtr) { return *(const U32*) memPtr; }
U32 MEM_readLE32(const void* memPtr)
{
	if (MEM_isLittleEndian())
		return MEM_read32(memPtr);
	else {
		const BYTE* p = (const BYTE*)memPtr;
		return (U32)((U32)p[0] + ((U32)p[1]<<8) + ((U32)p[2]<<16) + ((U32)p[3]<<24));
	}
}

int zstd_skipframe(ZSTDres *res) {
#ifdef DEBUG
debug_printprogress(res,"skipframe");
#endif
	if (res->inframepos == 0) {
		if (!zstd_readheader(res)) return(0);
	}
	if (res->zfh.frameType == ZSTD_skippableFrame) {
/*		size_t framesize = ZSTD_findFrameCompressedSize(res->buffer, res->inframepos);*/
		/* on update to zstd1.4.4, ZSTD_findFrameCompressedSize got framesize wrong
		 * just a quick hack to solve
		*/
		size_t framesize = 8 + MEM_readLE32((BYTE const*)res->buffer+4);
		DPRINT("skip %llu of skippable frame at ftello: %zu ----",framesize - res->inframepos, ftello(res->finput));
		fseeko(res->finput,framesize - res->inframepos,SEEK_CUR);
		res->inframepos = 0;
		res->frameread = 0;
	} else {
		DPRINT("skip frame at ftello: %zu ----\n",ftello(res->finput));
		res->bufferpos = res->inframepos;
		res->lastblock = 0;
		while (1) {
			uint32_t blocksize;
			blocksize = zstd_readblockheader(res);
			if (!blocksize) break;
			if (res->lastblock && res->zfh.checksumFlag) {
				blocksize += 4;
			}
			fseeko(res->finput,blocksize,SEEK_CUR);
			if (res->lastblock) break;
		}
		res->inframepos = 0;
		res->frameread = 0;
	}
	return(1);
}

ZSTDres *zstdopen(FILE* finput,FILE *findex) {
	ZSTDres *res = malloc(sizeof(ZSTDres));
	res->finput = finput;
	res->findex = findex;
	res->inframepos = 0;
	res->framepos = 0;
	res->contentsize = 0;
	res->usize = 0;
	res->indexbsize = 0;
	res->currentblock = BUFFER_EMPTY;
	res->currentpos = 0;
	res->buffer = (char *)malloc(ZSTD_FRAMEHEADERSIZE_MAX);
	res->buffersize = ZSTD_FRAMEHEADERSIZE_MAX;
	res->outbuffersize = 1024;
	res->outbuffer = (char *)malloc(res->outbuffersize);
	/* Read and Check Archive Header */
	zstd_readheader(res);

	if (findex != NULL) {
		char *buffer;
		size_t n = 40;
		ssize_t read;
		NODPRINT("indexfile used");
		/* open index file */
		res->findex = findex;
		buffer = (char *)malloc(40);
		read = getline(&buffer,&n,findex);
		if (read < 9 || strncmp(buffer,"#bym zsti",9) != 0) {
			fprintf(stderr,"not a zstindex file");
			exit(1);
		}
		fseeko(findex,0,SEEK_SET);
		while (1) {
			read = getline(&buffer,&n,findex);
			if (read == -1) {
				fprintf(stderr,"error in header of zstd index file\n");
				exit(1);
			}
			if (read == 4 && strncmp(buffer,"...",3) == 0) break;
			if (read < 6) continue;
			if (read > 7 && strncmp(buffer,"usize: ",7) == 0) {
				sscanf(buffer+7,"%llu",(long long unsigned int *)&(res->usize));
			} else if (read > 9 && strncmp(buffer,"binsize: ",9) == 0) {
				sscanf(buffer+9,"%llu",(long long unsigned int *)&(res->indexbsize));
			} else if (read > 11 && strncmp(buffer,"framesize: ",11) == 0) {
				sscanf(buffer+11,"%llu",(long long unsigned int *)&(res->framesize));
			} else if (read > 6 && strncmp(buffer,"name: ",6) == 0) {
				yaml_checkparam(buffer,read,6,"zstindex",8);
			} else if (read > 9 && strncmp(buffer,"version: ",9) == 0) {
				yaml_checkparam(buffer,read,9,"0.1",3);
			} else if (read > 10 && strncmp(buffer,"datatype: ",10) == 0) {
				yaml_checkparam(buffer,read,10,"uint64",6);
			} else if (read > 6 && strncmp(buffer,"type: ",6) == 0) {
				yaml_checkparam(buffer,read,6,"array",5);
			} else if (read > 11 && strncmp(buffer,"byteorder: ",11) == 0) {
				yaml_checkparam(buffer,read,11,"l",1);
			}
		}
		if (res->indexbsize == 0) {
			fprintf(stderr,"zstdi error: binsize missing or 0\n");
			exit(1);
		}
		if (res->usize == 0) {
			fprintf(stderr,"zstdi error: usize missing or 0\n");
			exit(1);
		}
		res->indexstart = ftell(findex);
		NODPRINT("indexfile usize: %llu",(long long unsigned int)res->usize);
		NODPRINT("indexfile indexstart: %llu",(long long unsigned int)res->indexstart);
		NODPRINT("indexfile indexbsize: %llu",(long long unsigned int)res->indexbsize);
	}
	return res;
}

void zstdclose(ZSTDres *res) {
	free(res->buffer);
	free(res->outbuffer);
	FCLOSE(res->finput);
	if (res->findex != NULL) FCLOSE(res->findex);
	free(res);
}

void zstdindex_write(FILE *findex,off_t pos) {
	putc_unlocked((unsigned char) pos,findex);
	putc_unlocked((unsigned char) (pos >> 8),findex);
	putc_unlocked((unsigned char) (pos >> 16),findex);
	putc_unlocked((unsigned char) (pos >> 24),findex);
	putc_unlocked((unsigned char) (pos >> 32),findex);
	putc_unlocked((unsigned char) (pos >> 40),findex);
	putc_unlocked((unsigned char) (pos >> 48),findex);
	putc_unlocked((unsigned char) (pos >> 56),findex);
}

uint64_t zstdindex_read(FILE *findex) {
	uint64_t pos;
	pos = (((uint64_t)getc_unlocked(findex)))
		| (((uint64_t)getc_unlocked(findex)) << 8)
		| (((uint64_t)getc_unlocked(findex)) << 16)
		| (((uint64_t)getc_unlocked(findex)) << 24)
		| (((uint64_t)getc_unlocked(findex)) << 32)
		| (((uint64_t)getc_unlocked(findex)) << 40)
		| (((uint64_t)getc_unlocked(findex)) << 48)
		| (((uint64_t)getc_unlocked(findex)) << 56);
	return(pos);
}

ZSTDres *zstd_openfile(char *file) {
	ZSTDres *res;
	FILE *finput, *findex;
	char *indexfile;
	finput = fopen64_or_die(file, "r");
	indexfile = zstd_findindex(file);
	findex = fopen(indexfile, "r");
	res = zstdopen(finput,findex);
	return(res);
}

int zstd_seek(ZSTDres *res, uint64_t pos, int where) {
	DPRINT("zstd_seek: %lu where:%d\n",pos,where);
	if (where == SEEK_SET) {
	} else if (where == SEEK_CUR) {
		pos = res->currentpos + pos;
	} else {
		fprintf(stderr,"internal error in zstd_seek\n");
		exit(1);
	}
	/* goto block containing currentpos */
	if (pos < res->framepos) {
		fseeko(res->finput,0,SEEK_SET);
		res->inframepos = 0;
		res->frameread = 0;
		res->contentsize = 0;
		res->currentpos = 0;
		res->framepos = 0;
		zstd_readheader(res);
	} else if (pos < (res->framepos + res->contentsize)) {
		res->currentpos = pos;
		return(1);
	}
	if (res->findex != NULL) {
		/* find block using the index */
		unsigned long long startblock, filepos;
		startblock = pos/res->framesize;
		NODPRINT("indexpos: %llu",(long long unsigned int)(res->indexstart + 8*startblock));
		fseeko(res->findex,res->indexstart + 8*startblock,SEEK_SET);
		res->inframepos = 0;
		res->frameread = 0;
		res->contentsize = 0;
		res->framepos = startblock * res->framesize;
		filepos = zstdindex_read(res->findex);
		NODPRINT("filepos: %llu",(long long unsigned int)filepos);
		fseeko(res->finput,filepos,SEEK_SET);
		zstd_readheader(res);
		res->currentpos = pos;
		return(1);
	} else {
		/* find block by skipping */
		if (res->inframepos == 0 && res->frameread == 1) {
			zstd_readheader(res);
		}
		while (1) {
			if (pos < res->framepos) {
				fprintf(stderr,"seek back not supported yet\n");
				exit(1);
			} else if (pos < (res->framepos + res->contentsize)) {
				res->currentpos = pos;
				return(1);
			}
#ifdef DEBUG
debug_printprogress(res,"seeking");
#endif
			zstd_skipframe(res);
			zstd_readheader(res);
			if (feof(res->finput)) break;
		}
		return(0);
	}
}

int zstd_read(ZSTDres *res, void *data, uint64_t size) {
	char *writebuffer, *dest = (char *)data;
	uint64_t read;
	unsigned int skip, writesize;
	skip = res->currentpos - res->framepos;
	if (res->contentsize == 0) {
		while (res->contentsize == 0) {
			if (feof(res->finput)) break;
			zstd_skipframe(res);
			zstd_readheader(res);
		}
		zstd_readframe(res);
	} else if (!res->frameread) {
		zstd_readframe(res);
	} else if (skip > res->contentsize) {
		zstd_readheader(res);
		while (res->contentsize == 0) {
			if (feof(res->finput)) break;
			zstd_skipframe(res);
			zstd_readheader(res);
		}
		zstd_readframe(res);
	}
	while (1) {
		if (size == 0) break;
		writebuffer = res->outbuffer + skip;
		writesize = res->contentsize - skip;
		skip = 0;
		if (writesize > size) {
			writesize = size; size = 0;
		} else {
			size -= writesize;
		}
		memcpy(dest,writebuffer,writesize);
		dest += writesize;
		if (size == 0) break;
		zstd_readheader(res);
		while (res->contentsize == 0) {
			if (feof(res->finput)) break;
			zstd_skipframe(res);
			zstd_readheader(res);
		}
		if (feof(res->finput)) break;
		zstd_readframe(res);
	}
	read = (dest - (char *)data);
	res->currentpos += read;
	return(read);
}

int zstd_get(ZSTDres *res) {
	char c;
	int read = zstd_read(res, &c, 1);
/* fprintf(stderr,"%c",c);fflush(stderr); */
	if (read == 0 && feof(res->finput)) {
		return(-1);
	} else {
		return((int)c);
	}
}
