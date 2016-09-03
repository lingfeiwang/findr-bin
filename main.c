/* Copyright 2016 Lingfei Wang
 * 
 * This file is part of Findr.
 * 
 * Findr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Findr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with Findr.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <base/logger.h>
#include <base/os.h>
#include <pij/gassist/gassist.h>
#include <pij/rank.h>
#include <base/random.h>
#include <base/macros.h>
#include <base/lib.h>

#define MACROSTR(X)	#X
#define STR(X)	MACROSTR(X)

#define	SLIBINFONAME	STR(LIBINFONAME)
#define	SLIBINFOVERSION	STR(LIBINFOVERSION)
#define	LOGLV_AUTO	6

//Checks binary and library versions are the same. Otherwise prints warning
void checkversion()
{
	const char bln[]=SLIBINFONAME;
	const char blv[]=SLIBINFOVERSION;
	const char* ln=lib_name();
	const char* lv=lib_version();
	char diff=0;
	
	if(strcmp(bln,ln))
	{
		diff=1;
		fprintf(stderr,"Warning: different names for binary interface and library.");
	}
	if(strcmp(blv,lv))
	{
		diff=1;
		fprintf(stderr,"Warning: different versions for binary interface and library.");
	}

	fprintf(stderr,"Binary interface for library %s, version %s.%s",bln,blv,_NEWLINE_);
	if(diff)
		fprintf(stderr,"Using library name %s, version %s.%s",ln,lv,_NEWLINE_);
}

//Display usage info.
void usage(FILE* fp,const char* binpath)
{
	fprintf(fp,"Usage: %s loglv rs nth method [method-args...]%s",binpath,_NEWLINE_);
	fprintf(fp,"\tloglv:\tLog level, 0 for automatic (warnings + errors)%s",_NEWLINE_);
	fprintf(fp,"\trs:\tInitial random seed. 0 for using current time.%s",_NEWLINE_);
	fprintf(fp,"\tnth:\tMaximum number of parallel threads, 0 for automatic%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"For each method, look up the function of the same name in the library headers for usage.%s",_NEWLINE_);
	fprintf(fp,"Supports the following methods:%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"pij_gassist_a: invokes function pijs_gassist_a in library%s",_NEWLINE_);
	fprintf(fp,"pij_gassist_tot: invokes function pijs_gassist_tot in library%s",_NEWLINE_);
	fprintf(fp,"Usage: pij_gassist_??? fg ft ft2 nt nt2 ns fp na nodiag%s",_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfg:\tPath of input file containing genotype data for E(A). Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp:\tPath of output file for probability of recommended test combinations. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tna:\tNumber of allleles for the species considered, = n_v-1.%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"pijs_gassist_a: invokes function pijs_gassist_a in library%s",_NEWLINE_);
	fprintf(fp,"pijs_gassist_tot: invokes function pijs_gassist_tot in library%s",_NEWLINE_);
	fprintf(fp,"Usage: pijs_gassist_??? fg ft ft2 nt nt2 ns fp1 fp2 fp3 fp4 fp5 na nodiag%s",_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfg:\tPath of input file containing genotype data for E(A). Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp1:\tPath of output file for probability of test 1, E->A v.s. E  A. Since significant eQTL inputs are expected, outputs in this file are manually set to all 1. Data vector size: (nt).%s",_NEWLINE_);
	fprintf(fp,"\tfp2:\tPath of output file for probability of test 2, E->B v.s. E  B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp3:\tPath of output file for probability of test 3, E->A->B v.s. E->A->B with E->B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp4:\tPath of output file for probability of test 4, E->A->B with E->B v.s. E->A  B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp5:\tPath of output file for probability of test 5, E->A->B with E->B v.s. A<-E->B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tna:\tNumber of allleles for the species considered, = n_v-1.%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"pij_rank_a: invokes function pij_rank_a in library%s",_NEWLINE_);
	fprintf(fp,"Usage: pij_rank_a ft ft2 nt nt2 ns fp nodiag%s",_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp:\tPath of output file for probability. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"File formats:%s",_NEWLINE_);
	fprintf(fp,"This binary interface accepts two input/output file formats: raw and tsv. Within each run, all input and output files must all have the same format, either raw or tsv. This is indicated from the method name. Normal method names indicate raw file format, and appending the method name with '_tsv' or '_csv' would indicate tsv file format.%s",_NEWLINE_);
	fprintf(fp,"\tRaw format: for matrices, row-major sequence is used. All genotype data follow type 'GTYPE', and all expression data and output probabilities follow type 'FTYPE'. Detailed definitions of GTYPE and FTYPE can be found in base/types.h in source tree, and make/def.dev in base folder. By default, GTYPE is 8-bit (1-byte) unsigned char, FTYPE is 32-bit (4-byte) signed float.%s",_NEWLINE_);
	fprintf(fp,"\tTsv format: tab separated value files. Each row should be separated by new line, and each column by tab or space. The number of rows and columns must match input dimensional parameters.%s",_NEWLINE_);
}

//vector and matrix tsv output functions
#define	OFORMAT	"%.8G"
static int vectorf_fprintf(FILE* fp,const VECTORF* v)
{
	size_t	i;
	if(!fprintf(fp,OFORMAT,VECTORFF(get)(v,0)))
		return 1;
	for(i=1;i<v->size;i++)
	{
		if(!putc(' ',fp))
			return 1;
		if(!fprintf(fp,OFORMAT,VECTORFF(get)(v,i)))
			return 1;
	}
	return !fputs(_NEWLINE_,fp);
}

static int matrixf_fprintf(FILE* fp,const MATRIXF* m)
{
	size_t	i;
	for(i=0;i<m->size1;i++)
	{
		VECTORFF(const_view)	vv=MATRIXFF(const_row)(m,i);
		if(vectorf_fprintf(fp,&vv.vector))
			return 1;
	}
	return 0;
}
#undef	OFORMAT

//Library initialization
void bin_call_lib_init(const char* argv[])
{
	unsigned char loglv;
	unsigned long rs;
	size_t nth;
	loglv=(unsigned char)atoi(argv[0]);
	sscanf(argv[1],"%lu",&rs);
	nth=(size_t)atoi(argv[2]);
	if(!loglv)
		loglv=LOGLV_AUTO;
	lib_init(loglv,rs,nth);
}


/*****************************************************************************
 *	Exposed functions
 ****************************************************************************/

//pij_rank generic function, for raw and tsv versions
int bin_pij_rank_func(int argc,const char* argv[],int (*func)(const MATRIXF*,const MATRIXF*,MATRIXF*,char),const char funcname[],int (*fin_fm)(FILE*,MATRIXF*),int (*fout_fm)(FILE*,const MATRIXF*))
{
#define CLEANUP CLEANMATF(t)CLEANMATF(t2)CLEANMATF(p)\
				if(fp){fclose(fp);fp=0;}
	const char *f_t,*f_t2,*f_p;
	char nodiag;
	size_t	ns,ng,nt;
	FILE	*fp=0;
	MATRIXF	*t=0,*t2=0,*p=0;
	
	if(argc!=7)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_t=argv[0];
	f_t2=argv[1];
	ng=(size_t)atoi(argv[2]);
	nt=(size_t)atoi(argv[3]);
	ns=(size_t)atoi(argv[4]);
	if(!(ng&&nt&&ns))
	{
		LOG(0,"Invalid input dimensions.")
		return -1;
	}
	f_p=argv[5];
	nodiag=(char)atoi(argv[6]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value.")
		return -1;
	}
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	t=MATRIXFF(alloc)(ng,ns);
	t2=MATRIXFF(alloc)(nt,ns);
	p=MATRIXFF(alloc)(ng,nt);
	if(!(t&&t2&&p))
		ERRRET("Not enough memory.")
	
	//File reads
	LOG(11,"Reading file %s.",f_t)
	fp=fopen(f_t,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_t)
	if(fin_fm(fp,t))
		ERRRET("Can't read file or has wrong size: %s.",f_t)
	fclose(fp);	
	LOG(11,"Reading file %s.",f_t2)
	fp=fopen(f_t2,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_t2)
	if(fin_fm(fp,t2))
		ERRRET("Can't read file or has wrong size: %s.",f_t2)
	fclose(fp);
	fp=0;
	
	//Calculation
	if(func(t,t2,p,nodiag))
		ERRRET("%s failed.",funcname)
	
	//File writes
	LOG(11,"Writing file %s.",f_p)
	fp=fopen(f_p,"wb");
	if(!fp)
		ERRRET("Can't open file %s.", f_p)
	if(fout_fm(fp,p))
		ERRRET("Can't write to file %s.",f_p)
	fclose(fp);
	fp=0;
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_pij_rank_a(int argc,const char* argv[])
{
	return bin_pij_rank_func(argc,argv,pij_rank_a,"pij_rank_a",MATRIXFF(fread),MATRIXFF(fwrite));
}

static inline int bin_pij_rank_a_tsv(int argc,const char* argv[])
{
	return bin_pij_rank_func(argc,argv,pij_rank_a,"pij_rank_a_tsv",MATRIXFF(fscanf),matrixf_fprintf);
}

//pijs_gassist generic function, for pijs_gassist_a and pijs_gassist_tot, for raw and tsv versions
int bin_pijs_gassist_func(int argc,const char* argv[],int (*func)(const MATRIXG*,const MATRIXF*,const MATRIXF*,VECTORF*,MATRIXF*,MATRIXF*,MATRIXF*,MATRIXF*,size_t,char),const char funcname[],int (*fin_gm)(FILE*,MATRIXG*),int (*fin_fm)(FILE*,MATRIXF*),int (*fout_fv)(FILE*,const VECTORF*),int (*fout_fm)(FILE*,const MATRIXF*))
{
#define CLEANUP CLEANMATG(g)CLEANMATF(t)CLEANMATF(t2)CLEANVECF(p1)\
				CLEANMATF(p2)CLEANMATF(p3)CLEANMATF(p4)CLEANMATF(p5)\
				if(fp){fclose(fp);fp=0;}
#define FWRITE(F,FP,P)	LOG(11,"Writing file %s.",FP)\
						fp=fopen(FP,"wb");\
						if(!fp)\
							ERRRET("Can't open file %s.", FP)\
						if(F(fp,P))\
							ERRRET("Can't write to file %s.",FP)\
						fclose(fp);
	const char *f_g,*f_t,*f_t2,*f_p1,*f_p2,*f_p3,*f_p4,*f_p5;
	char nodiag;
	size_t	ns,ng,nt;
	size_t	nv;
	FILE	*fp=0;
	VECTORF		*p1=0;
	MATRIXG		*g=0;
	MATRIXF		*t=0,*t2=0,*p2=0,*p3=0,*p4=0,*p5=0;
	
	if(argc!=13)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[0];
	f_t=argv[1];
	f_t2=argv[2];
	ng=(size_t)atoi(argv[3]);
	nt=(size_t)atoi(argv[4]);
	ns=(size_t)atoi(argv[5]);
	f_p1=argv[6];
	f_p2=argv[7];
	f_p3=argv[8];
	f_p4=argv[9];
	f_p5=argv[10];
	nv=(size_t)atoi(argv[11])+1;
	if(!(ng&&nt&&ns&&(nv-1)))
	{
		LOG(0,"Invalid input dimensions or allele count.")
		return -1;
	}
	nodiag=(char)atoi(argv[12]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value.")
		return -1;
	}
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	g=MATRIXGF(alloc)(ng,ns);
	t=MATRIXFF(alloc)(ng,ns);
	t2=MATRIXFF(alloc)(nt,ns);
	p1=VECTORFF(alloc)(ng);
	p2=MATRIXFF(alloc)(ng,nt);
	p3=MATRIXFF(alloc)(ng,nt);
	p4=MATRIXFF(alloc)(ng,nt);
	p5=MATRIXFF(alloc)(ng,nt);
	if(!(g&&t&&t2&&p1&&p2&&p3&&p4&&p5))
		ERRRET("Not enough memory.")
	
	//File reads
	LOG(11,"Reading file %s.",f_g)
	fp=fopen(f_g,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_g)
	if(fin_gm(fp,g))
		ERRRET("Can't read file or has wrong size: %s.",f_g)
	fclose(fp);
	LOG(11,"Reading file %s.",f_t)
	fp=fopen(f_t,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_t)
	if(fin_fm(fp,t))
		ERRRET("Can't read file or has wrong size: %s.",f_t)
	fclose(fp);	
	LOG(11,"Reading file %s.",f_t2)
	fp=fopen(f_t2,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_t2)
	if(fin_fm(fp,t2))
		ERRRET("Can't read file or has wrong size: %s.",f_t2)
	fclose(fp);
	fp=0;
	
	//Calculation
	if(func(g,t,t2,p1,p2,p3,p4,p5,nv,nodiag))
		ERRRET("%s failed.",funcname)
	
	//File writes
	FWRITE(fout_fv,f_p1,p1)
	FWRITE(fout_fm,f_p2,p2)
	FWRITE(fout_fm,f_p3,p3)
	FWRITE(fout_fm,f_p4,p4)
	FWRITE(fout_fm,f_p5,p5)
	fp=0;
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef FWRITE
#undef CLEANUP
}

//pij_gassist generic function, for pij_gassist_a and pij_gassist_tot, for raw and tsv versions
int bin_pij_gassist_func(int argc,const char* argv[],int (*func)(const MATRIXG*,const MATRIXF*,const MATRIXF*,MATRIXF*,size_t,char),const char funcname[],int (*fin_gm)(FILE*,MATRIXG*),int (*fin_fm)(FILE*,MATRIXF*),int (*fout_fm)(FILE*,const MATRIXF*))
{
#define CLEANUP CLEANMATG(g)CLEANMATF(t)CLEANMATF(t2)CLEANMATF(p)\
				if(fp){fclose(fp);fp=0;}
#define FWRITE(F,FP,P)	LOG(11,"Writing file %s.",FP)\
						fp=fopen(FP,"wb");\
						if(!fp)\
							ERRRET("Can't open file %s.", FP)\
						if(F(fp,P))\
							ERRRET("Can't write to file %s.",FP)\
						fclose(fp);
	const char *f_g,*f_t,*f_t2,*f_p;
	char nodiag;
	size_t	ns,ng,nt;
	size_t	nv;
	FILE	*fp=0;
	MATRIXG		*g=0;
	MATRIXF		*t=0,*t2=0,*p=0;
	
	if(argc!=9)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[0];
	f_t=argv[1];
	f_t2=argv[2];
	ng=(size_t)atoi(argv[3]);
	nt=(size_t)atoi(argv[4]);
	ns=(size_t)atoi(argv[5]);
	f_p=argv[6];
	nv=(size_t)atoi(argv[7])+1;
	if(!(ng&&nt&&ns&&(nv-1)))
	{
		LOG(0,"Invalid input dimensions or allele count.")
		return -1;
	}
	nodiag=(char)atoi(argv[8]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value.")
		return -1;
	}
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	g=MATRIXGF(alloc)(ng,ns);
	t=MATRIXFF(alloc)(ng,ns);
	t2=MATRIXFF(alloc)(nt,ns);
	p=MATRIXFF(alloc)(ng,nt);
	if(!(g&&t&&t2&&p))
		ERRRET("Not enough memory.")
	
	//File reads
	LOG(11,"Reading file %s.",f_g)
	fp=fopen(f_g,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_g)
	if(fin_gm(fp,g))
		ERRRET("Can't read file or has wrong size: %s.",f_g)
	fclose(fp);
	LOG(11,"Reading file %s.",f_t)
	fp=fopen(f_t,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_t)
	if(fin_fm(fp,t))
		ERRRET("Can't read file or has wrong size: %s.",f_t)
	fclose(fp);	
	LOG(11,"Reading file %s.",f_t2)
	fp=fopen(f_t2,"rb");
	if(!fp)
		ERRRET("Can't open file %s.",f_t2)
	if(fin_fm(fp,t2))
		ERRRET("Can't read file or has wrong size: %s.",f_t2)
	fclose(fp);
	fp=0;
	
	//Calculation
	if(func(g,t,t2,p,nv,nodiag))
		ERRRET("%s failed.",funcname)
	
	//File writes
	FWRITE(fout_fm,f_p,p)
	fp=0;
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef FWRITE
#undef CLEANUP
}

static inline int bin_pijs_gassist_a(int argc,const char* argv[])
{
	return bin_pijs_gassist_func(argc,argv,pijs_gassist_a,"pijs_gassist_a",MATRIXGF(fread),MATRIXFF(fread),VECTORFF(fwrite),MATRIXFF(fwrite));
}

static inline int bin_pijs_gassist_tot(int argc,const char* argv[])
{
	return bin_pijs_gassist_func(argc,argv,pijs_gassist_tot,"pijs_gassist_tot",MATRIXGF(fread),MATRIXFF(fread),VECTORFF(fwrite),MATRIXFF(fwrite));
}

static inline int bin_pijs_gassist_a_tsv(int argc,const char* argv[])
{
	return bin_pijs_gassist_func(argc,argv,pijs_gassist_a,"pijs_gassist_a_tsv",MATRIXGF(fscanf),MATRIXFF(fscanf),vectorf_fprintf,matrixf_fprintf);
}

static inline int bin_pijs_gassist_tot_tsv(int argc,const char* argv[])
{
	return bin_pijs_gassist_func(argc,argv,pijs_gassist_tot,"pijs_gassist_tot_tsv",MATRIXGF(fscanf),MATRIXFF(fscanf),vectorf_fprintf,matrixf_fprintf);
}

static inline int bin_pij_gassist_a(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist_a,"pij_gassist_a",MATRIXGF(fread),MATRIXFF(fread),MATRIXFF(fwrite));
}

static inline int bin_pij_gassist_tot(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist_tot,"pij_gassist_tot",MATRIXGF(fread),MATRIXFF(fread),MATRIXFF(fwrite));
}

static inline int bin_pij_gassist_a_tsv(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist_a,"pij_gassist_a_tsv",MATRIXGF(fscanf),MATRIXFF(fscanf),matrixf_fprintf);
}

static inline int bin_pij_gassist_tot_tsv(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist_tot,"pij_gassist_tot_tsv",MATRIXGF(fscanf),MATRIXFF(fscanf),matrixf_fprintf);
}

int main(int argc,const char* argv[])
{
#define	GETFUNCNAMERAW(X)	if(!strcmp(argv[4],STR(X))) func=bin_##X;
#define	GETFUNCNAMETSV(X)	if(!strcmp(argv[4],STR(X)"_tsv")) func=bin_##X##_tsv;
#define	GETFUNCNAMECSV(X)	if(!strcmp(argv[4],STR(X)"_csv")) func=bin_##X##_tsv;
#define	GETFUNCNAMES(X)		GETFUNCNAMERAW(X) else GETFUNCNAMETSV(X) \
							else GETFUNCNAMECSV(X)
	int	ret;
	int	(*func)(int,const char*[]);
	
	checkversion();
	func=0;
	//Look for method name
	if(argc>=5)
	{
		GETFUNCNAMES(pijs_gassist_a)
		else GETFUNCNAMES(pijs_gassist_tot)
		else GETFUNCNAMES(pij_gassist_a)
		else GETFUNCNAMES(pij_gassist_tot)
		else GETFUNCNAMES(pij_rank_a)
	}
	//If method name not found
	if(!func)
	{
		usage(stdout,argv[0]);
		return 0;
	}
	bin_call_lib_init(argv+1);
	ret=func(argc-5,argv+5);
	//If wrong number of arguments
	return ret;
}
#undef	GETFUNCNAME
#undef	LIBINFO


































