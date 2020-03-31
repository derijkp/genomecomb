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
#include <inttypes.h>
#include "tools.h"
#include "debug.h"


/*
 *----------------------------------------------------------------------
 *
 * NeedReversing --
 *
 *  Taken mostly from Tcl code of the same name
 *
 * Results:
 *	0	No re-ordering needed.
 *	1	Reverse the bytes:	01234567 <-> 76543210 (little to big)
 *	2	Apply this re-ordering: 01234567 <-> 45670123 (Nokia to little)
 *	3	Apply this re-ordering: 01234567 <-> 32107654 (Nokia to big)
 *
 * Side effects:
 *	None
 *
 *----------------------------------------------------------------------
 */
static int NeedReversing(int format) {
    switch (format) {
	/* native floats and doubles: never reverse */
    case 'd':
    case 'f':
	/* big endian ints: never reverse */
    case 'I':
    case 'S':
    case 'W':
#ifdef WORDS_BIGENDIAN
	/* native ints: reverse if we're little-endian */
    case 'n':
    case 't':
    case 'm':
	/* f: reverse if we're little-endian */
    case 'Q':
    case 'R':
#else /* !WORDS_BIGENDIAN */
	/* small endian floats: reverse if we're big-endian */
    case 'r':
#endif /* WORDS_BIGENDIAN */
	return 0;

#ifdef WORDS_BIGENDIAN
	/* small endian floats: reverse if we're big-endian */
    case 'q':
    case 'r':
#else /* !WORDS_BIGENDIAN */
	/* native ints: reverse if we're little-endian */
    case 'n':
    case 't':
    case 'm':
	/* f: reverse if we're little-endian */
    case 'R':
#endif /* WORDS_BIGENDIAN */
	/* small endian ints: always reverse */
    case 'i':
    case 's':
    case 'w':
	return 1;
    }
    return 0;
}

/*
 *----------------------------------------------------------------------
 *
 * CopyNumber --
 * also taken from Tcl
 *
 *	This routine is called by FormatNumber and ScanNumber to copy a
 *	floating-point number. If required, bytes are reversed while copying.
 *	The behaviour is only fully defined when used with IEEE float and
 *	double values (guaranteed to be 4 and 8 bytes long, respectively.)
 *
 * Results:
 *	None
 *
 * Side effects:
 *	Copies length bytes
 *
 *----------------------------------------------------------------------
 */

static void
CopyNumber(
    const void *from,		/* source */
    void *to,			/* destination */
    unsigned int length,	/* Number of bytes to copy */
    int type)			/* What type of thing are we copying? */
{
    switch (NeedReversing(type)) {
    case 0: 
	memcpy(to, from, length);
	break;
    case 1: {
	const unsigned char *fromPtr = from;
	unsigned char *toPtr = to;

	switch (length) {
	case 4:
	    toPtr[0] = fromPtr[3];
	    toPtr[1] = fromPtr[2];
	    toPtr[2] = fromPtr[1];
	    toPtr[3] = fromPtr[0];
	    break;
	case 8:
	    toPtr[0] = fromPtr[7];
	    toPtr[1] = fromPtr[6];
	    toPtr[2] = fromPtr[5];
	    toPtr[3] = fromPtr[4];
	    toPtr[4] = fromPtr[3];
	    toPtr[5] = fromPtr[2];
	    toPtr[6] = fromPtr[1];
	    toPtr[7] = fromPtr[0];
	    break;
	}
	break;
    }
    case 2: {
	const unsigned char *fromPtr = from;
	unsigned char *toPtr = to;

	toPtr[0] = fromPtr[4];
	toPtr[1] = fromPtr[5];
	toPtr[2] = fromPtr[6];
	toPtr[3] = fromPtr[7];
	toPtr[4] = fromPtr[0];
	toPtr[5] = fromPtr[1];
	toPtr[6] = fromPtr[2];
	toPtr[7] = fromPtr[3];
	break;
    }
    case 3: {
	const unsigned char *fromPtr = from;
	unsigned char *toPtr = to;

	toPtr[0] = fromPtr[3];
	toPtr[1] = fromPtr[2];
	toPtr[2] = fromPtr[1];
	toPtr[3] = fromPtr[0];
	toPtr[4] = fromPtr[7];
	toPtr[5] = fromPtr[6];
	toPtr[6] = fromPtr[5];
	toPtr[7] = fromPtr[4];
	break;
    }
    }
}

int bcol_read(char *ftype,long double *result) {
	int buffer[9],type;
	long long value;
	float fvalue;
	double dvalue;
	int i, v;
	type = ftype[0];
	switch(type) {
	case 'c':
		v = fgetc(stdin);
		if (v == EOF) {return 0;}
		buffer[0] = v;
		value = buffer[0];
		if ((ftype[1] != 'u')) {
			if (value & 0x80) {
				value |= -0x100;
			}
		}
		*result = (long double)value;
		break;
	case 's':
		for(i = 0 ; i < 2 ; i++) {
			v = fgetc(stdin);
NODPRINT("%d",v)
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		if (NeedReversing(type)) {
			value = (long long) (buffer[0] + (buffer[1] << 8));
		} else {
			value = (long long) (buffer[1] + (buffer[0] << 8));
		}
NODPRINT("%d,%d -> %" PRId64 "",buffer[0],buffer[1],value)
		if ((ftype[1] != 'u')) {
			if (value & 0x8000) {
				value |= -0x10000;
			}
		}
NODPRINT("   -> %" PRId64 "",value)
		*result = (long double)value;
		break;
	case 'i':
		for(i = 0 ; i < 4 ; i++) {
			v = fgetc(stdin);
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		if (NeedReversing(type)) {
		value = (long long) (buffer[0]
			+ (buffer[1] << 8)
			+ (buffer[2] << 16)
			+ (((long long)buffer[3]) << 24));
		} else {
		value = (long long) (buffer[3]
			+ (buffer[2] << 8)
			+ (buffer[1] << 16)
			+ (((long long)buffer[0]) << 24));
		}
NODPRINT("%d,%d,%d,%d -> %" PRId64 "",buffer[0],buffer[1],buffer[2],buffer[3],value)
		if ((ftype[1] != 'u')) {
			if ((value & (((unsigned int)1)<<31)) && (value > 0)) {
			    value -= (((unsigned int)1)<<31);
			    value -= (((unsigned int)1)<<31);
			}
		}
		*result = (long double)value;
		break;
	case 'w':
		for(i = 0 ; i < 8 ; i++) {
			v = fgetc(stdin);
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		if (NeedReversing(type)) {
		value = ((long long) buffer[0])
			| (((long long) buffer[1]) << 8)
			| (((long long) buffer[2]) << 16)
			| (((long long) buffer[3]) << 24)
			| (((long long) buffer[4]) << 32)
			| (((long long) buffer[5]) << 40)
			| (((long long) buffer[6]) << 48)
			| (((long long) buffer[7]) << 56);
		} else {
		value = ((long long) buffer[7])
			| (((long long) buffer[6]) << 8)
			| (((long long) buffer[5]) << 16)
			| (((long long) buffer[4]) << 24)
			| (((long long) buffer[3]) << 32)
			| (((long long) buffer[2]) << 40)
			| (((long long) buffer[1]) << 48)
			| (((long long) buffer[0]) << 56);
		}
		*result = (long double)value;
		break;
	case 'f':
		for(i = 0 ; i < 4 ; i++) {
			v = fgetc(stdin);
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		CopyNumber(buffer, &fvalue, sizeof(float), type);
		*result = (long double)fvalue;
		break;
	case 'd':
		for(i = 0 ; i < 8 ; i++) {
			v = fgetc(stdin);
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		CopyNumber(buffer, &dvalue, sizeof(double), type);
		*result = (long double)dvalue;
		break;
	}
	return 1;
}

#define DBL_MAX      1.79769313486231470e+308

int main(int argc, char *argv[]) {
	char *format,*chr;
	long double value,min,max;
	int shift,pos,accept;
	int status,begin=0,start,error;
	if ((argc != 7)) {
		fprintf(stderr,"Format is: getregions chromosome format start min max shift");
		exit(EXIT_FAILURE);
	}
	chr = argv[1];
	format = argv[2];
	start = atoi(argv[3]);
	min = atoll(argv[4]);
	if (argv[5][0] == '\0') {
		max = DBL_MAX;
	} else {
		max = atoll(argv[5]);
	}
	shift = atoi(argv[6]);
	pos = start - shift;
	begin = -1;
	status = 0;
	while (1) {
		error = bcol_read(format,&value);
NODPRINT("%d -> %Lf",pos,value)
		if (error == 0) break;
		accept = (value >= min && value <= max);
		if (status) {
			if (!accept) {
				fprintf(stdout,"%s\t%d\t%d\n", chr, begin+shift, pos+shift);
				status = 0;
			}
		} else {
			if (accept) {
				begin = pos;
				status = 1;
			}
		}
		pos++;
	}
	if (status) {
		fprintf(stdout,"%s\t%d\t%d\n", chr, begin+shift, pos+shift);
	}
	exit(EXIT_SUCCESS);
}
