/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"
#include "debug.h"
#include "hash.h"

DString *extractID(DString *string,DString *id) {
	char *str, *end;
	str = strstr(string->string,"ID=");
	if (str == NULL) {
		DStringSetS(id,"",0);
		return id;
	}
	str+= 3;
	end = strstr(str,",");
	DStringSetS(id,str,end-str);
	return id;
}

int main(int argc, char *argv[]) {
	Hash_table *conv_formata;
	FILE *fd = NULL,*fo = NULL;
	DString *line=NULL, *string=NULL, *temp=NULL;
	DString *snp=DStringNewFromChar("snp"), *del=DStringNewFromChar("del"), *ins=DStringNewFromChar("ins"), *sub=DStringNewFromChar("sub");
	DString *geno = NULL, *outinfo = NULL;
	DStringArray *header=NULL, *format=NULL, *info=NULL, *samples=NULL;
	DStringArray *formatfields=NULL, *headerfields=NULL, *infofields=NULL, *linea=NULL;
	DString *ref=DStringNew(), *alt=DStringNew(), *id=DStringNew();
	int *order = NULL;
	int read,i,j,pos,maxtab,igeno,isample,linenr=0;
	line = DStringNew();
	if (argc >= 2) {fd = fopen(argv[1],"r");} else {fd = stdin;}
	if (argc == 3) {fo = fopen(argv[2],"w");} else {fo = stdout;}
	NODPRINT("==== Reading header ====")
	DStringGetLine(line, fd);
	if (line->size < 16 || strncmp(line->string,"##fileformat=VCF",16) != 0) {
		fprintf(stderr,"error: input is not a vcf file\n");
		exit(1);
	}
	temp = DStringNew();
	conv_formata = hash_init();
	dstring_hash_set(conv_formata,DStringNewFromChar("AD"),(void *)DStringNewFromChar("alleledepth"));
	dstring_hash_set(conv_formata,DStringNewFromChar("GT"),(void *)DStringNewFromChar("genotype"));
	dstring_hash_set(conv_formata,DStringNewFromChar("DP"),(void *)DStringNewFromChar("coverage"));
	dstring_hash_set(conv_formata,DStringNewFromChar("FT"),(void *)DStringNewFromChar("filter"));
	dstring_hash_set(conv_formata,DStringNewFromChar("GL"),(void *)DStringNewFromChar("loglikelihood"));
	dstring_hash_set(conv_formata,DStringNewFromChar("GQ"),(void *)DStringNewFromChar("genoqual"));
	dstring_hash_set(conv_formata,DStringNewFromChar("HQ"),(void *)DStringNewFromChar("haploqual"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AN"),(void *)DStringNewFromChar("totalallelecount"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AC"),(void *)DStringNewFromChar("allelecount"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AF"),(void *)DStringNewFromChar("frequency"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AA"),(void *)DStringNewFromChar("Ancestralallele"));
	dstring_hash_set(conv_formata,DStringNewFromChar("DB"),(void *)DStringNewFromChar("dbsnp"));
	dstring_hash_set(conv_formata,DStringNewFromChar("H2"),(void *)DStringNewFromChar("Hapmap2"));
	fprintf(fo,"# -- tsv converted from vcf, original comments follow --\n");
	info = DStringArrayNew(10);
	format = DStringArrayNew(10);
	while ((read = DStringGetLine(line, fd)) != -1) {
		linenr++;
		fprintf(fo,"%s\n",line->string);
		if (line->string[0] == '#') {
			if (line->string[1] == '#') {
				if (line->size > 9 && strncmp("FORMAT=",line->string+2,7) == 0) {
					DStringArrayAppend(format,line->string+9,-1);
				} else if (line->size > 7 && strncmp("INFO=",line->string+2,5) == 0) {
					DStringArrayAppend(info,line->string+7,-1);
				}
			} else {
				header = DStringArrayFromCharM(line->string+1," \t");
				break;
			}
		}
	}
	fprintf(fo,"# ----\n");
	fprintf(fo,"chromosome\tbegin\tend\ttype\tref\talt\tname\tquality\tfilter");
	if (header->size >= 10) {
		samples = DStringArrayRange(header,9,header->size-1);
	} else {
		samples = DStringArrayNew(0);
	}
	if (header->size >= 9) {
		NODPRINT("==== Parsing format ====")
		formatfields = DStringArrayNew(10);
		DStringArrayAppend(formatfields,"GT",2);
		headerfields = DStringArrayNew(10);
		DStringArrayAppend(headerfields,"alleleSeq1",10);
		DStringArrayAppend(headerfields,"alleleSeq2",10);
		DStringArrayAppend(headerfields,"zyg",3);
		DStringArrayAppend(headerfields,"phased",6);
		DStringArrayAppend(headerfields,"genotypes",9);
		for (i = 0 ; i < format->size ; i ++) {
			DString *ds;
			id = extractID(DStringArrayGet(format,i),id);
			if (strcmp(id->string,"GT") == 0) continue;
			DStringArrayAppend(formatfields,id->string,id->size);
			ds = (DString *)dstring_hash_get(conv_formata,id);
			if (ds != NULL) {
				DStringArrayAppend(headerfields,ds->string,ds->size);
			} else {
				DStringArrayAppend(headerfields,id->string,id->size);
			}
		}
		NODPRINT("==== Parsing header/samples ====")
		if (samples->size <= 1) {
			for (i = 0 ; i < headerfields->size ; i++) {
				DString *temp = DStringArrayGet(headerfields,i);
				fprintf(fo,"\t%*.*s",temp->size,temp->size,temp->string);
			}
		} else {
			for (j = 0 ; j < samples->size ; j++) {
				for (i = 0 ; i < headerfields->size ; i++) {
					DString *ds = DStringArrayGet(headerfields,i);
					DStringSetS(temp, ds->string, ds->size);
					DStringAppendS(temp, "-", 1);
					DStringAppendS(temp, DStringArrayGet(samples,j)->string,DStringArrayGet(samples,j)->size);
					fprintf(fo,"\t%*.*s",temp->size,temp->size,temp->string);
				}
			}
		}
	}
	NODPRINT("\n\n==== Parsing info ====")
	infofields = DStringArrayNew(10);
	for (i = 0 ; i < info->size ; i ++) {
		DString *ds;
		id = extractID(DStringArrayGet(info,i),id);
		DStringArrayAppend(infofields,id->string,id->size);
		if (id->string[0] == 'D' && id->string[1] == 'P' && id->string[2] == '\0') {
			fprintf(fo,"\ttotalcoverage");
		} else {
			ds = (DString *)dstring_hash_get(conv_formata,id);
			if (ds != NULL) {
				fprintf(fo,"\t%*.*s",ds->size,ds->size,ds->string);
			} else {
				fprintf(fo,"\t%*.*s",id->size,id->size,id->string);
			}
		}
	}
	fprintf(fo,"\n");
	maxtab = 9+samples->size;
#define a_chrom(array) (DStringArrayGet(array,0))
#define a_pos(array) (DStringArrayGet(array,1))
#define a_id(array) (DStringArrayGet(array,2))
#define a_ref(array) (DStringArrayGet(array,3))
#define a_alt(array) (DStringArrayGet(array,4))
#define a_qual(array) (DStringArrayGet(array,5))
#define a_filter(array) (DStringArrayGet(array,6))
#define a_info(array) (DStringArrayGet(array,7))
#define a_format(array) (DStringArrayGet(array,8))
#define a_geno(array) (DStringArrayGet(array,9))
	outinfo =  (DString *)malloc(infofields->size*sizeof(DString));
	linea = DStringArrayNew(maxtab+2);
	NODPRINT("==== Parsing data ====")
	while ((read = DStringGetTab(line,fd,maxtab,linea,1,NULL)) != -1) {
		DStringArray *lineformat,*alts;
		DString *type;
		char *genotypecur,zyg;
		int l1,l2,begin,end,len,diff;
		linenr++;
		lineformat = DStringArrayFromChar(a_format(linea)->string,':');
		/* set genos [lrange $line 9 end] */
		l1 = a_ref(linea)->size;
		l2 = 0;
		alts = DStringArrayFromChar(a_alt(linea)->string,',');
		for (i = 0 ; i < alts->size ; i ++) {
			DString *temp = DStringArrayGet(alts,i);
			if (temp->size > l2) {l2 = temp->size;}
		}
		pos = atoi(a_pos(linea)->string);
		if (l1 == 1 && l2 == 1) {
			diff = 0;
			begin = pos - 1;
			end = pos;
			type = snp;
		} else {
			diff = 1;
			begin = pos;
			end = pos + l1 - 1;
			if (l1 == 1) {
				type = ins;
			} else if (l2 == 1) {
				type = del;
			} else {
				type = sub;
			}
		}
		if (l1 > 20) {
			DStringPrintf(ref,"%d",l1 - 1);
		} else {
			DStringSetS(ref,a_ref(linea)->string+diff,a_ref(linea)->size-diff);
		}
		DStringSetS(alt,DStringArrayGet(alts,0)->string+diff,DStringArrayGet(alts,0)->size-diff);
		for (i = 1; i < alts->size ; i++) {
			DString *temp = DStringArrayGet(alts,i);
			DStringAppendS(alt, ",",1);
			DStringAppendS(alt, temp->string+diff, temp->size-diff);
		}
		fprintf(fo,"%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s", a_chrom(linea)->string, begin, end, type->string, ref->string, alt->string, a_id(linea)->string, a_qual(linea)->string, a_filter(linea)->string);
		if (header->size >= 9) {
			NODPRINT("==== Process genos ====")
			order = realloc(order,formatfields->size*sizeof(int));
			for (i = 0 ; i < formatfields->size ; i++) {
				order[i] = DStringArraySearch(lineformat,DStringArrayGet(formatfields,i)->string,DStringArrayGet(formatfields,i)->size);
			}
			igeno = 9; /* genos start at col 9 */
			for (isample = 0 ; isample < samples->size ; isample++) {
				DStringArray *genoa;
				DString *genotype, *temp;
				int phased,i,a1,a2;
				geno = DStringArrayGet(linea,igeno++);
				genoa = DStringArrayFromChar(geno->string,':');
				if (order[0] == -1) {
					genotype = DStringEmtpy();
				} else {
					genotype = DStringArrayGet(genoa,order[0]);
				}
				genotypecur = genotype->string;
				while (*genotypecur != '|' && *genotypecur != '/' && *genotypecur != '\0') {
					genotypecur++;
				}
				if (*genotypecur == '/') {phased = 0;} else {phased = 1;}
				if (genotype->string[0] == '.') {
					fprintf(fo,"\t?");
					a1 = -2;
				} else if (genotype->string == genotypecur) {
					fprintf(fo,"\t-");
					a1 = -1;
				} else {
					a1 = atol(genotype->string);
					if (a1 == 0) {
						fprintf(fo,"\t%*.*s",ref->size,ref->size,ref->string);
					} else {
						temp = DStringArrayGet(alts,a1-1);
						fprintf(fo,"\t%*.*s",temp->size-diff,temp->size-diff,temp->string+diff);
					}
				}
				if (*genotypecur != '\0') {genotypecur++;}
				if (*genotypecur == '.') {
					fprintf(fo,"\t?");
					a2 = -2;
				} else if (*genotypecur == '\0') {
					fprintf(fo,"\t-");
					a2 = -1;
				} else {
					a2 = atol(genotypecur);
					if (a2 == 0) {
						fprintf(fo,"\t%*.*s",ref->size,ref->size,ref->string);
					} else {
						temp = DStringArrayGet(alts,a2-1);
						fprintf(fo,"\t%*.*s",temp->size-diff,temp->size-diff,temp->string+diff);
					}
				}
				if (a1 != 0) {
					if (a2 == a1) {
						zyg = 'm';
					} else if (a2 != 0) {
						zyg = 'c';
					} else {
						zyg = 't';
					}
				} else if (a2 != 0) {
					zyg = 't';
				} else {
					zyg = 'r';
				}
				fprintf(fo,"\t%c",zyg);
				fprintf(fo,"\t%d",phased);
				genotypecur = genotype->string;
				/* fprintf(fo,"\t%s",genotype->string); */
				fputc_unlocked('\t',fo);
				while (*genotypecur != '\0') {
					if (*genotypecur == '|') {
						fputc_unlocked(',',fo);
					} else if (*genotypecur == '/') {
						fputc_unlocked(';',fo);
					} else {
						fputc_unlocked(*genotypecur,fo);
					}
					genotypecur++;
				}
				for (i = 1 ; i < formatfields->size; i++) {
					if (order[i] <= 0) {
						fprintf(fo,"\t");
					} else {
						fprintf(fo,"\t%s",DStringArrayGet(genoa,order[i])->string);
					}
				}
				if (genoa != NULL) DStringArrayDestroy(genoa);
			}
		}
		/* info output */
		{
			char *cur,*curend;
			DString *info = a_info(linea);
			for (i = 0 ; i< infofields->size ; i++) {
				outinfo[i].size = -1;
			}
			cur = info->string;
			curend = cur;
			while (1) {
				while (*curend != '=' && *curend != ';' && *curend != '\0') curend++;
				if (curend == cur) {
					cur = ++curend;
					continue;
				}
				pos = DStringArraySearch(infofields,cur,curend-cur);
				if (pos == -1) {fprintf(stderr,"line %d: info field %*.*s not described in header, skipping\n",linenr,(int)(curend-cur),(int)(curend-cur),cur);}
				if (*curend == '=') {curend++;}
				cur = curend;
				if (*curend) {while (*curend != ';' && *curend != '\0') curend++;}
				outinfo[pos].string = cur;
				outinfo[pos].size = curend-cur;
				if (*curend == '\0') break;
				cur = ++curend;
			}
			for (i = 0 ; i< infofields->size ; i++) {
				len  = outinfo[i].size;
				if (len == -1) {
					fprintf(fo,"\t");
				} else if (len == 0) {
					fprintf(fo,"\t1");
				} else {
					fprintf(fo,"\t%*.*s",len,len,outinfo[i].string);
				}
			}
			fprintf(fo,"\n");
		}
		DStringArrayDestroy(lineformat);
		DStringArrayDestroy(alts);
	}
	if (fd != stdin) {fclose(fd);}
	/* free dstrings */
	DStringDestroy(line);
	DStringDestroy(snp);
	DStringDestroy(ins);
	DStringDestroy(del);
	DStringDestroy(sub);
	DStringDestroy(string);
	DStringDestroy(temp);
	DStringDestroy(ref);
	DStringDestroy(alt);
	DStringDestroy(id);
	/* free arrays */
	DStringArrayDestroy(linea);
	DStringArrayDestroy(format);
	DStringArrayDestroy(samples);
	DStringArrayDestroy(formatfields);
	DStringArrayDestroy(headerfields);
	hash_destroy(conv_formata,(hash_free_func *)DStringDestroy,(hash_free_func *)DStringDestroy);
	/* free other */
	free(order);
	exit(EXIT_SUCCESS);
}
