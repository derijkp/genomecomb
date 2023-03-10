#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <inttypes.h>
#include "tools.h"
#include "lz4tools.h"

#ifdef __MINGW32__
#include "getline.c"
#endif

unsigned LZ4IO_readLE32 (const void* s) {
    const unsigned char* srcPtr = (const unsigned char*)s;
    unsigned value32 = srcPtr[0];
    value32 += (srcPtr[1]<<8);
    value32 += (srcPtr[2]<<16);
    value32 += ((unsigned)srcPtr[3])<<24;
    return value32;
}

char *lz4_findindex(char *filename) {
	char *outfile = (char *)malloc(strlen(filename) + 6);
	sprintf(outfile,"%s.lz4i",filename);
	return outfile;
}

#define yaml_checkparam(buffer,read,start,test,size) { \
if (read < (start+size) || strncmp(buffer+start,test,size) != 0 || (buffer[start+size] != ' ' && buffer[start+size] != '\n')) { \
	fprintf(stderr,"lz4i error: only %s supported\n",test); \
	exit(1); \
} \
}

LZ4res *lz4open(FILE* finput,FILE *findex,char *file) {
	LZ4res *res = malloc(sizeof(LZ4res));
	unsigned char MNstore[MAGICNUMBER_SIZE];
	unsigned magicNumber;
	BYTE FLG, BD, HC;
	uint64_t SIZE;
	unsigned blockSizeID;

	res->finput = finput;
	res->findex = findex;
	res->writebuffer = NULL;
	res->writesize = 0;
	/* Check Archive Header */
	READ(finput,MAGICNUMBER_SIZE,&MNstore,"empty file")
	magicNumber = LZ4IO_readLE32(MNstore);   /* Little Endian format */
	if (magicNumber == LEGACY_MAGICNUMBER) {
		if (file != NULL) fprintf(stderr,"file %s ",file);
		ERROR("unsupported format, legacy lz4 format not supported");
	} else if (magicNumber == LZ4IO_SKIPPABLE0) {
		if (file != NULL) fprintf(stderr,"file %s ",file);
		ERROR("unsupported format, start with skippable frame not supported");
	} else if (magicNumber != LZ4IO_MAGICNUMBER) {
		if (file != NULL) fprintf(stderr,"file %s ",file);
		ERROR("unsupported format, not lz4");
	}
	res->usize = 0;
	res->indexbsize = 0;
	res->currentblock = BUFFER_EMPTY;
	res->currentpos = 0;
	/* read flags */
	READ(finput,1,&FLG,"lz4 header error")
	res->version = (FLG>>6) & _2BITS;
	res->blockMode = (FLG>>5) & _1BIT;
	if (res->blockMode != 1) ERROR("Files with dependend blocks are not supported")
	res->blockChecksumFlag = (FLG>>4) & _1BIT;
	res->contentSizeFlag = (FLG>>3) & _1BIT;
	res->contentChecksumFlag = (FLG>>2) & _1BIT;
	
	/* read blockSizeID */
	READ(finput,1,&BD,"lz4 header error")
	blockSizeID = (BD>>4) & _3BITS;
	switch (blockSizeID) {
		case 4: res->blocksize = 65536; break;
		case 5: res->blocksize = 262144; break;
		case 6: res->blocksize = 1048576; break;
		case 7: res->blocksize = 4194304; break;
		default: ERROR("unsupported blocksize");
	}
	res->buffer = (char *)malloc(res->blocksize+4);
	res->outbuffer = (char *)malloc(res->blocksize);
	/* content size present? */
	if (res->contentSizeFlag) {
		READ(finput,8,&SIZE,"lz4 header error")
	}

	/* Header Checksum */
	READ(finput,1,&HC,"lz4 header error")
	res->startblocks = (int)ftello(finput);

	/* print info */
	NODPRINT("version: %d",res->version);
	NODPRINT("blockMode: %d",res->blockMode);
	NODPRINT("blockChecksumFlag: %d",res->blockChecksumFlag);
	NODPRINT("blockSizeID: %d",blockSizeID);
	NODPRINT("contentChecksumFlag: %d",res->contentChecksumFlag);

	if (findex != NULL) {
		char *buffer;
		size_t n = 40;
		ssize_t read;
		NODPRINT("indexfile used");
		/* open index file */
		res->findex = findex;
		buffer = (char *)malloc(40);
		read = getline(&buffer,&n,findex);
		if (read < 9 || strncmp(buffer,"#bym lz4i",9) != 0) {
			fprintf(stderr,"not a lz4index file");
			exit(1);
		}
		fseeko(findex,0,SEEK_SET);
		while (1) {
			read = getline(&buffer,&n,findex);
			if (read == -1) {
				fprintf(stderr,"error in header of lz4 index file\n");
				exit(1);
			}
			if (read == 4 && strncmp(buffer,"...",3) == 0) break;
			if (read < 6) continue;
			if (read > 7 && strncmp(buffer,"usize: ",7) == 0) {
				sscanf(buffer+7,"%" PRIu64 "",(uint64_t *)&(res->usize));
			} else if (read > 9 && strncmp(buffer,"binsize: ",9) == 0) {
				sscanf(buffer+9,"%" PRIu64 "",(uint64_t *)&(res->indexbsize));
			} else if (read > 6 && strncmp(buffer,"name: ",6) == 0) {
				yaml_checkparam(buffer,read,6,"lz4index",8);
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
			fprintf(stderr,"lz4i error: binsize missing or 0\n");
			exit(1);
		}
		if (res->usize == 0) {
			fprintf(stderr,"lz4i error: usize missing or 0\n");
			exit(1);
		}
		res->indexstart = ftell(findex);
		NODPRINT("indexfile usize: %" PRIu64 "",res->usize);
		NODPRINT("indexfile indexstart: %" PRIu64 "",res->indexstart);
		NODPRINT("indexfile indexbsize: %" PRIu64 "",res->indexbsize);
	}
	return res;
}

void lz4close(LZ4res *res) {
	free(res->buffer);
	free(res->outbuffer);
	FCLOSE(res->finput);
	if (res->findex != NULL) FCLOSE(res->findex);
	free(res);
}

int lz4_readblock(LZ4res *res, unsigned int startblock) {
	FILE *finput = res->finput, *findex = res->findex;
	uint64_t pos;
	unsigned int compressedsize, count;
	int readsize = 0, numdecompressed, uncompressed, r;
	if (startblock == res->currentblock) {
		return 1;
	}
	NODPRINT("res->currentblock: %d",res->currentblock);
	NODPRINT("startblock: %d",startblock);
	if (res->writebuffer == NULL || startblock != res->currentblock+1) {
		if (findex != NULL) {
			/* find block using the index */
			NODPRINT("indexpos: %" PRIu64 "",(res->indexstart + 8*startblock));
			fseeko(findex,res->indexstart + 8*startblock,SEEK_SET);
			pos = lz4index_read(findex);
			NODPRINT("filepos: %" PRIu64 "",pos);
			fseeko(finput,pos,SEEK_SET);
		} else {
			/* find block by skipping */
			count = startblock;
			if (res->writebuffer == NULL || startblock < res->currentblock) {
				fseeko(finput,res->startblocks,SEEK_SET);
			} else {
				count = count - res->currentblock - 1;
			}
			while (1) {
				if (!count) break;
				NODPRINT("pos:       %ld",ftello(finput));
				READ(finput,4,res->buffer,"could not read block size")
				compressedsize = LZ4IO_readLE32(res->buffer);
				if (compressedsize == 0) break;
				uncompressed = compressedsize >= 2147483648u;
				if (uncompressed) {
					compressedsize -= 2147483648u;
				}
				NODPRINT("compressedsize: %d (%d)",compressedsize,uncompressed);
				if (res->blockChecksumFlag) {
					readsize = compressedsize+4;
				} else {
					readsize = compressedsize;
				}
				fseeko(finput,readsize,SEEK_CUR);
				count--;
			}
		}
	}
	NODPRINT("pos:       %ld",ftello(finput));
	READ(finput,4,res->buffer,"could not read block size")
	compressedsize = LZ4IO_readLE32(res->buffer);
	if (compressedsize == 0) {return 0;}
	uncompressed = compressedsize >= 2147483648u;
	if (uncompressed) {
		compressedsize -= 2147483648u;
	}
	NODPRINT("compressedsize: %d (%d)",compressedsize,uncompressed);
	if (res->blockChecksumFlag) {
		readsize = compressedsize+4;
	} else {
		readsize = compressedsize;
	}
	r = fread(res->buffer,1,readsize,finput);
	if (r != readsize) {
		if (r < 0) {ERROR("error reading lz4 file");}
		ERROR("lz4 file prematurely ends");
	}
	if (uncompressed) {
		res->writebuffer = res->buffer;
		res->writesize = compressedsize;
	} else {
		numdecompressed = LZ4_decompress_safe (res->buffer, res->outbuffer, compressedsize, res->blocksize);
		NODPRINT("numdecompressed: %d",numdecompressed);
		res->writebuffer = res->outbuffer;
		res->writesize = numdecompressed;
	}
	res->currentblock = startblock;
	NODPRINT("after read res->currentblock: %d",res->currentblock);
	return 1;
}

void lz4index_write(FILE *findex,off_t pos) {
	putc_unlocked((unsigned char) pos,findex);
	putc_unlocked((unsigned char) (pos >> 8),findex);
	putc_unlocked((unsigned char) (pos >> 16),findex);
	putc_unlocked((unsigned char) (pos >> 24),findex);
	putc_unlocked((unsigned char) (pos >> 32),findex);
	putc_unlocked((unsigned char) (pos >> 40),findex);
	putc_unlocked((unsigned char) (pos >> 48),findex);
	putc_unlocked((unsigned char) (pos >> 56),findex);
}

uint64_t lz4index_read(FILE *findex) {
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

LZ4res *lz4_openfile(char *file, int useindex) {
	LZ4res *res;
	FILE *finput, *findex;
	char *indexfile;
	finput = fopen64_or_die(file, "r");
	if (useindex) {
		indexfile = lz4_findindex(file);
		findex = fopen64_or_die(indexfile, "r");
	} else {
		findex = NULL;
	}
	res = lz4open(finput,findex,file);
	return(res);
}

void lz4_seek(LZ4res *res, int64_t pos, int where) {
	if (where == SEEK_SET) {
		res->currentpos = pos;
	} else if (where == SEEK_CUR) {
		if (pos < 0 && (uint64_t)(-pos) > res->currentpos) {
			res->currentpos = 0;
		} else {
			res->currentpos += pos;
		}
	} else {
		fprintf(stderr,"internal error in lz4_seek\n");
		exit(1);
	}
}

int lz4_read(LZ4res *res, void *data, uint64_t size) {
	char *writebuffer, *dest = (char *)data;
	uint64_t startblock, read;
	unsigned int skip, writesize;
	startblock = res->currentpos/res->blocksize;
	skip = res->currentpos-(startblock*res->blocksize);
	/* decompress */
	while (1) {
		if (size == 0) break;
		if (!lz4_readblock(res,startblock++)) break;
		writebuffer = res->writebuffer + skip;
		writesize = res->writesize - skip;
		skip = 0;
		if (writesize > size) {
			writesize = size; size = 0;
		} else {
			size -= writesize;
		}
		memcpy(dest,writebuffer,writesize);
		dest += writesize;
	}
	read = (dest - (char *)data);
	res->currentpos += read;
	return(read);
}

int lz4_get(LZ4res *res) {
	uint64_t startblock;
	unsigned int skip;
	startblock = res->currentpos/res->blocksize;
	skip = res->currentpos-(startblock*res->blocksize);
	/* decompress */
	if (!lz4_readblock(res,startblock++) || skip >= res->writesize) {return EOF;}
	res->currentpos++;
	return (res->writebuffer[skip]);
}
