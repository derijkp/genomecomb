#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "razf.h"

#define WINDOW_SIZE 4096

static int razf_main_usage()
{
	printf("\n");
	printf("Usage:   razip [options] [file] ...\n\n");
	printf("Options: -c      write on standard output, keep original files unchanged\n");
	printf("         -d      decompress\n");
	printf("         -l      list compressed file contents\n");
	printf("         -b INT  decompress at INT position in the uncompressed file\n");
	printf("         -s INT  decompress INT bytes in the uncompressed file\n");
	printf("         -h      give this help\n");
	printf("\n");
	return 0;
}

static int write_open(const char *fn, int is_forced)
{
	int fd = -1;
	char c;
	if (!is_forced) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC | O_EXCL, 0644)) < 0 && errno == EEXIST) {
			printf("razip: %s already exists; do you wish to overwrite (y or n)? ", fn);
			scanf("%c", &c);
			if (c != 'Y' && c != 'y') {
				printf("razip: not overwritten\n");
				exit(1);
			}
		}
	}
	if (fd < 0) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0) {
			fprintf(stderr, "razip: %s: Fail to write\n", fn);
			exit(1);
		}
	}
	return fd;
}

void razip_compress (int f_src, int f_dst) {
	int c;
	RAZF *rz;
	void *buffer;
	int64_t start, end, size;
	start = 0; size = -1; end = -1;
	rz = razf_dopen(f_dst, "w");
	buffer = malloc(WINDOW_SIZE);
	while((c = read(f_src, buffer, WINDOW_SIZE)) > 0) razf_write(rz, buffer, c);
	razf_close(rz); // f_dst will be closed here
	free(buffer);
	close(f_src);
}

int razip_decompress (char *filename, char *dst, int64_t start, int64_t end, int is_forced) {
	int c, f_dst;
	RAZF *rz;
	char *name;
	void *buffer;
	if (dst != NULL) {
		char *ext = dst+strlen(dst) - 3;
		if ((strncmp(ext,".rz",3) != 0) && (strncmp(ext,".gz",3) != 0)) {
			fprintf(stderr,"razip: %s: unknown suffix \"%s\" -- ignored\n", dst, ext);
			return 1;
		}
		name = strdup(dst);
		name[strlen(name) - 3] = '\0';
		f_dst = write_open(name, is_forced);
		free(name);
	} else {
		f_dst = fileno(stdout);
	}
	rz = razf_open(filename, "r");
	buffer = malloc(WINDOW_SIZE);
	razf_seek(rz, start, SEEK_SET);
	while(1){
		if(end < 0) c = razf_read(rz, buffer, WINDOW_SIZE);
		else c = razf_read(rz, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
		if(c <= 0) break;
		start += c;
		write(f_dst, buffer, c);
		if(end >= 0 && start >= end) break;
	}
	razf_close(rz);
	free(buffer);
	return 0;
}

int main(int argc, char **argv)
{
	int c, compress, pstdout, is_forced, error;
	RAZF *rz;
	void *buffer;
	int64_t start, end, size;

	compress = 1; pstdout = 0; start = 0; size = -1; end = -1; is_forced = 0;
	while((c  = getopt(argc, argv, "cdlhfb:s:")) >= 0){
		switch(c){
		case 'h': return razf_main_usage();
		case 'd': compress = 0; break;
		case 'c': pstdout = 1; break;
		case 'l': compress = 2; break;
		case 'b': start = atoll(optarg); break;
		case 's': size = atoll(optarg); break;
		case 'f': is_forced = 1; break;
		}
	}
	if (size >= 0) end = start + size;
	if(end >= 0 && end < start){
		fprintf(stderr, " -- Illegal region: [%ld, %ld] --\n", start, end);
		return 1;
	}
	if(compress == 1){
		int f_src, f_dst = -1;
		if(argc > optind){
			while (optind < argc) {
				char *ext = argv[optind]+strlen(argv[optind]) - 3;
				if ((strncmp(ext,".rz",3) == 0) || (strncmp(ext,".gz",3) == 0)) {
					fprintf(stderr,"razip: %s: already has suffix \"%s\" for compressed data -- ignored\n", argv[optind], ext);
					optind++; continue;
				}
				if((f_src = open(argv[optind], O_RDONLY)) < 0){
					fprintf(stderr, " -- Cannot open file: %s --\n", argv[optind]);
					return 1;
				}
				if(pstdout){
					f_dst = fileno(stdout);
				} else {
					char *name = NULL;
					name = malloc((strlen(argv[optind]) + 5)*sizeof(char));
					strcpy(name, argv[optind]);
					strcat(name, ".rz");
					f_dst = write_open(name, is_forced);
					if (f_dst < 0) goto error;
					free(name);
				}
				razip_compress(f_src,f_dst);
				if (!pstdout) unlink(argv[optind]);
				close(f_src);
				optind++;
			}
			return 0;
		} else if (pstdout) { 
			f_src = fileno(stdin);
			f_dst = fileno(stdout);
			razip_compress(f_src,f_dst);
			close(f_src);
			return 0;
		} else {
			return razf_main_usage();
		}
	} else {
		if(argc <= optind) return razf_main_usage();
		if(compress == 2){
			rz = razf_open(argv[optind], "r");
			if(rz->file_type == FILE_TYPE_RZ) {
				printf("%20s%20s%7s %s\n", "compressed", "uncompressed", "ratio", "name");
				printf("%20lld%20lld%6.1f%% %s\n", (long long)rz->end, (long long)rz->src_end, rz->end * 100.0f / rz->src_end,argv[optind]);
			} else {
				fprintf(stdout, "%s is not a regular rz file\n", argv[optind]);
			}
		} else {
			int f_dst;
			if (!pstdout) {
				char *name;
/*
				if (strstr(argv[optind], ".rz") - argv[optind] != strlen(argv[optind]) - 3) {
					printf("razip: %s: unknown suffix -- ignored\n", argv[optind]);
					goto error;
				}
*/
				while (optind < argc) {
					error = razip_decompress(argv[optind],argv[optind],start,end,is_forced);
					if (!error && !pstdout) unlink(argv[optind]);
					optind++;
				}
			} else {
				razip_decompress(argv[optind],NULL,start,end,is_forced);
			}
		}
		return 0;
	}
	error:
		return 1;

}
