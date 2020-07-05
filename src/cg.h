/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#define _GNU_SOURCE

#include <stdio.h>
#include <stdint.h>

#ifdef __MINGW32__
#define ssize_t long int
#define fopen64(filename,mode) (fopen(filename,mode))
#define putc_unlocked(c, stream) (putc(c, stream))
#define getc_unlocked(stream) (getc(stream))
#define fputc_unlocked(c, stream) (fputc(c, stream))
ssize_t getline(char **lineptr, size_t *n, FILE *stream);
#endif

#define FCLOSE(f) ({if (fclose(f) == EOF) {perror("Error writing file: "); exit(EXIT_FAILURE);}})
