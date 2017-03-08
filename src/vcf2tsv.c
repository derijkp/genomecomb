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

char extractNumber(DString *string) {
	char *str;
	str = strstr(string->string,"Number=");
	if (str == NULL) {
		return '.';
	}
	str+= 7;
	return str[0];
}

void getfield(DString *result,char *list, int pos) {
	char *cur = list;
	int size = 0;
	while (pos > 0) {
		if (*cur == '\0') break;
		if (*cur == ';') break;
		if (*cur == ',') {
			pos--;
		}
		cur++;
	}
	result->string = cur;
	while (*cur != '\0' && *cur != ',' && *cur != ';') {
		cur++;
		size++;
	}
	result->size = size;
}

char numberfromid(DString *id) {
	if (id->size >= 2 && (strncmp(id->string,"AD",2) == 0 || strncmp(id->string,"RP",2) == 0)) {
		return('R');
	}
	if (id->size == 2 && (strncmp(id->string,"AC",2) == 0 || strncmp(id->string,"AF",2) == 0)) {
		return('A');
	}
	return('.');
}

typedef struct altvar {
	DString *type;
	int begin;
	int end;
	char *ref;
	int refsize;
	char *alt;
	int altsize;
} AltVar;

int main(int argc, char *argv[]) {
	Hash_table *conv_formata;
	FILE *fd = NULL,*fo = NULL;
	AltVar *altvars = NULL;
	DString *line=NULL, *string=NULL, *temp=NULL;
	DString *snp=DStringNewFromChar("snp"), *del=DStringNewFromChar("del"), *ins=DStringNewFromChar("ins"), *sub=DStringNewFromChar("sub");
	DString *geno = NULL, *outinfo = NULL, *formatfieldsnumber=NULL, *infofieldsnumber=NULL;
	DStringArray *header=NULL, *format=NULL, *info=NULL, *samples=NULL;
	DStringArray *formatfields=NULL , *headerfields=NULL, *infofields=NULL, *linea=NULL;
	DString *ref=DStringNew(), *alt=DStringNew(), *id=DStringNew();
	DString *num=DStringNew();
	int *order = NULL;
	int read,i,j,pos,maxtab,igeno,isample,linenr=0,numalleles=1,curallele, altvarsmax=0;
	DStringSetS(num,".",1);
	altvars = (AltVar *)malloc(5*sizeof(AltVar)); altvarsmax = 5;
	line = DStringNew();
	if (argc < 1) {
		fprintf(stderr,"Format is: vcf2tsv ?infile? ?outfile?\n");
		exit(EXIT_FAILURE);
	}
	if (argc >= 2) {
		fd = fopen(argv[1],"r");
		if (fd == NULL) {fprintf(stderr,"file %s does not exists\n",argv[1]);exit(1);}
	} else {fd = stdin;}
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
	fprintf(fo,"#filetype\ttsv/varfile\n");
	fprintf(fo,"#fileversion\t%s\n",FILEVERSION);
	fprintf(fo,"#split\t1\n");
	fprintf(fo,"#info\ttsv converted from vcf, original comments follow\n");
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
		formatfieldsnumber = DStringNew();
		DStringAppendS(formatfieldsnumber,"G",1);
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
			num->string[0] = extractNumber(DStringArrayGet(format,i));
			if (num->string[0] == '.') {
				num->string[0] = numberfromid(id);
			}
			DStringAppendS(formatfieldsnumber,num->string,1);
			ds = (DString *)dstring_hash_get(conv_formata,id);
			if (ds == NULL) {ds = id;}
			if (num->string[0] == 'R') {
				DStringSetS(temp,ds->string,ds->size);
				DStringAppendS(temp,"_ref",4);
				DStringArrayAppend(headerfields,temp->string,temp->size);
			}
			DStringArrayAppend(headerfields,ds->string,ds->size);
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
	infofieldsnumber = DStringNew();
	for (i = 0 ; i < info->size ; i ++) {
		DString *ds;
		id = extractID(DStringArrayGet(info,i),id);
		DStringArrayAppend(infofields,id->string,id->size);
		num->string[0] = extractNumber(DStringArrayGet(info,i));
		if (num->string[0] == '.') {
			num->string[0] = numberfromid(id);
		}
		DStringAppendS(infofieldsnumber,num->string,1);
		if (id->string[0] == 'D' && id->string[1] == 'P' && id->string[2] == '\0') {
			fprintf(fo,"\ttotalcoverage");
		} else {
			ds = (DString *)dstring_hash_get(conv_formata,id);
			if (ds == NULL) {ds = id;}
			if (num->string[0] == 'R') {
				fprintf(fo,"\t%*.*s_ref",ds->size,ds->size,ds->string);
			}			
			fprintf(fo,"\t%*.*s",ds->size,ds->size,ds->string);
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
		char zyg;
		int l1,l2,len;
		linenr++;
		lineformat = DStringArrayFromChar(a_format(linea)->string,':');
		/* set genos [lrange $line 9 end] */
		alts = DStringArrayFromChar(a_alt(linea)->string,',');
		numalleles = alts->size;
		if (numalleles > altvarsmax) {
			altvars = (AltVar *)realloc(altvars,numalleles*sizeof(AltVar));
			altvarsmax = numalleles;
		}
		/* determine type, ref, alt, ... for different alleles */
		for (curallele = 0 ; curallele < numalleles; curallele++) {
			DString *altallele = DStringArrayGet(alts,curallele);
			AltVar *altvar = altvars+curallele;
			char *curref, *curalt;
			/* determine type for this altallele */
			pos = atoi(a_pos(linea)->string);
			curref = a_ref(linea)->string;
			l1 = a_ref(linea)->size;
			curalt = altallele->string;
			l2 = altallele->size;
			pos--; /* vcf 1 based */
			if (*curref == *curalt) {
				/* remove base before simple indel */
				pos++;
				curref++; l1--;
				curalt++; l2--;
			}
			/* remove extra bases at end (from other overlapping alleles) */
			while (l1 > 0 && l2 > 0 && curref[l1-1] == curalt[l2-1]) {
				l1--; l2--;
			}
			altvar->ref = curref; altvar->refsize = l1;
			altvar->alt = curalt; altvar->altsize = l2;
			if (l1 == 1 && l2 == 1) {
				altvar->type = snp;
				altvar->begin = pos;
				altvar->end = pos + 1;
			} else if (l1 == 0 && l2 != 0) {
				altvar->type = ins;
				altvar->begin = pos;
				altvar->end = pos;
			} else if (l1 != 0 && l2 == 0) {
				altvar->type = del;
				altvar->begin = pos;
				altvar->end = pos+l1;
			} else {
				altvar->type = sub;
				altvar->begin = pos;
				altvar->end = pos+l1;
			}
		}
		/* go over the different alleles */
		for (curallele = 1 ; curallele <= numalleles; curallele++) {
			AltVar *altvar = altvars+curallele-1;
			NODPRINT("==== Print variant info ====")
			fprintf(fo,"%s\t%d\t%d\t%s", a_chrom(linea)->string, altvar->begin, altvar->end, altvar->type->string);
			if (altvar->refsize >= 20) {
				fprintf(fo,"\t%d", altvar->refsize);
			} else {
				fprintf(fo,"\t%*.*s", altvar->refsize, altvar->refsize, altvar->ref);
			}
			fprintf(fo,"\t%*.*s", altvar->altsize, altvar->altsize, altvar->alt);
			fprintf(fo,"\t%s\t%s\t%s", a_id(linea)->string, a_qual(linea)->string, a_filter(linea)->string);
			if (header->size >= 9) {
				NODPRINT("==== Determine geno fields order ====")
				order = realloc(order,formatfields->size*sizeof(int));
				for (i = 0 ; i < formatfields->size ; i++) {
					order[i] = DStringArraySearch(lineformat,DStringArrayGet(formatfields,i)->string,DStringArrayGet(formatfields,i)->size);
				}
				NODPRINT("==== Process genos ====")
				igeno = 9; /* genos start at col 9 */
				for (isample = 0 ; isample < samples->size ; isample++) {
					DStringArray *genoa;
					DString *genotype, *genotypelist;
					char *genotypestring,*genotypecur;
					int phased,i,a1,a2;
					genotypelist = DStringNew();
					geno = DStringArrayGet(linea,igeno++);
					genoa = DStringArrayFromChar(geno->string,':');
					if (order[0] == -1) {
						genotypestring = "0/0";
					} else {
						genotype = DStringArrayGet(genoa,order[0]);
						genotypestring = genotype->string;
					}
					genotypecur = genotypestring;
					while (*genotypecur != '|' && *genotypecur != '/' && *genotypecur != '\0') {
						genotypecur++;
					}
					if (*genotypecur == '/') {phased = 0;} else {phased = 1;}
					if (genotypestring[0] == '.') {
						fprintf(fo,"\t?");
						a1 = -2;
						DStringAppendS(genotypelist,"?",1);
					} else if (genotypestring == genotypecur) {
						/* empty genotype */
						fprintf(fo,"\t-");
						a1 = -1;
						DStringAppendS(genotypelist,"?",1);
					} else {
						a1 = atol(genotypestring);
						if (a1 == 0) {
							if (altvar->refsize >= 20) {
								fprintf(fo,"\t%d", altvar->refsize);
							} else {
								fprintf(fo,"\t%*.*s", altvar->refsize, altvar->refsize, altvar->ref);
							}
							DStringAppendS(genotypelist,"0",1);
						} else if (a1 == curallele) {
							fprintf(fo,"\t%*.*s",altvar->altsize,altvar->altsize,altvar->alt);
							DStringAppendS(genotypelist,"1",1);
						} else {
							AltVar *temp = altvars+a1-1;
							if (temp->type == altvar->type && temp->begin == altvar->begin && temp->end == altvar->end) {
								fprintf(fo,"\t%*.*s",temp->altsize,temp->altsize,temp->alt);
							} else {
								fprintf(fo,"\t@");
							}
							DStringAppendS(genotypelist,"2",1);
						}
					}
					if (*genotypecur == '|') {
						DStringAppendS(genotypelist,",",1);
					} else if (*genotypecur == '/') {
						DStringAppendS(genotypelist,";",1);
					} else {
						DStringAppendS(genotypelist,"?",1);
					}
					if (*genotypecur != '\0') {genotypecur++;}
					if (*genotypecur == '.') {
						fprintf(fo,"\t?");
						a2 = -2;
						DStringAppendS(genotypelist,"?",1);
					} else if (*genotypecur == '\0') {
						fprintf(fo,"\t-");
						a2 = -1;
						DStringAppendS(genotypelist,"?",1);
					} else {
						a2 = atol(genotypecur);
						if (a2 == 0) {
							if (altvar->refsize >= 20) {
								fprintf(fo,"\t%d", altvar->refsize);
							} else {
								fprintf(fo,"\t%*.*s", altvar->refsize, altvar->refsize, altvar->ref);
							}
							DStringAppendS(genotypelist,"0",1);
						} else if (a2 == curallele) {
							fprintf(fo,"\t%*.*s",altvar->altsize,altvar->altsize,altvar->alt);
							DStringAppendS(genotypelist,"1",1);
						} else {
							AltVar *temp = altvars+a2-1;
							if (temp->type == altvar->type && temp->begin == altvar->begin && temp->end == altvar->end) {
								fprintf(fo,"\t%*.*s",temp->altsize,temp->altsize,temp->alt);
							} else {
								fprintf(fo,"\t@");
							}
							DStringAppendS(genotypelist,"2",1);
						}
					}
					if (a1 < 0 || a2 < 0) {
						zyg = '?';
					} else if (a1 == curallele) {
						if (a2 == a1) {
							zyg = 'm';
						} else if (a2 > 0) {
							zyg = 'c';
						} else {
							zyg = 't';
						}
					} else if (a2 == curallele) {
						if (a1 > 0) {
							zyg = 'c';
						} else {
							zyg = 't';
						}
					} else if (a1 > 0 || a2 > 0) {
						zyg = 'o';
					} else {
						zyg = 'r';
					}
					fprintf(fo,"\t%c",zyg);
					fprintf(fo,"\t%d",phased);
					fputc_unlocked('\t',fo);
					DStringputs(genotypelist,fo);
					while (*genotypecur >= 48 && *genotypecur <= 57) {
						genotypecur++;
					}
					while (*genotypecur != '\0') {
						if (*genotypecur == '|') {
							fputc_unlocked(',',fo);
							genotypecur++;
						} else if (*genotypecur == '/') {
							fputc_unlocked(';',fo);
							genotypecur++;
						} else {
							a1 = atol(genotypecur);
							if (a1 == 0) {
								fprintf(fo,"0");
							} else if (a1 == curallele) {
								fprintf(fo,"1");
							} else {
								fprintf(fo,"2");
							}
							while (*genotypecur >= 48 && *genotypecur <= 57) {
								genotypecur++;
							}
							fputc_unlocked(*genotypecur,fo);
						}
					}
					for (i = 1 ; i < formatfields->size; i++) {
						if (order[i] <= 0) {
							fprintf(fo,"\t");
							if (formatfieldsnumber->string[i] == 'R') {
								fprintf(fo,"\t");
							}
						} else if (formatfieldsnumber->string[i] == 'R') {
							DString result, *value;
							value = DStringArrayGet(genoa,order[i]);
							getfield(&result,value->string,0);
							fprintf(fo,"\t%*.*s",result.size,result.size,result.string);
							getfield(&result,value->string,curallele);
							fprintf(fo,"\t%*.*s",result.size,result.size,result.string);
						} else if (formatfieldsnumber->string[i] == 'A') {
							DString result, *value;
							value = DStringArrayGet(genoa,order[i]);
							getfield(&result,value->string,curallele-1);
							fprintf(fo,"\t%*.*s",result.size,result.size,result.string);
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
						if (infofieldsnumber->string[i] == 'R') {
							fprintf(fo,"\t");
						}
					} else if (len == 0) {
						if (infofieldsnumber->string[i] == 'R') {
							fprintf(fo,"\t");
						}
						fprintf(fo,"\t1");
					} else if (infofieldsnumber->string[i] == 'R') {
						DString result;
						getfield(&result,outinfo[i].string,0);
						fprintf(fo,"\t%*.*s",result.size,result.size,result.string);
						getfield(&result,outinfo[i].string,curallele);
						fprintf(fo,"\t%*.*s",result.size,result.size,result.string);
					} else if (infofieldsnumber->string[i] == 'A') {
						DString result;
						getfield(&result,outinfo[i].string,curallele-1);
						fprintf(fo,"\t%*.*s",result.size,result.size,result.string);
					} else {
						fprintf(fo,"\t%*.*s",len,len,outinfo[i].string);
					}
				}
				fprintf(fo,"\n");
			}
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
