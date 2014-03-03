typedef struct VariantPos {
	int chr;
	int start;
	int end;
	int type;
	int ref;
	int alt;
	int a1;
	int a2;
	int max;
} VariantPos;

typedef struct Variant {
	DString *chr;
	int start;
	int end;
	DString *type;
	DString *ref;
	DString *alt;
} Variant;

/* returns 0 if equal, 1 if equal not including alt, 2 if different */
int varchecksort(Variant *prev,Variant *var,char *filename,int *nextpos);
void varputs(Variant var,FILE *f);
void result2var(DStringArray *result,VariantPos varpos, Variant *var);
int varcompare(Variant *var1, Variant *var2, int split);
void varpos_init(VariantPos *varpos);
int varpos_max(VariantPos *varpos);
