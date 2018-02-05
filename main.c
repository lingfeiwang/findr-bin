/* Copyright 2016-2018 Lingfei Wang
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
#define	ENABLE_PV
#define	ENABLE_CASSIST
#define	ENABLE_NETR

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <base/logger.h>
#include <base/os.h>
#include <pij/gassist/gassist.h>
#ifdef	ENABLE_CASSIST
#include <pij/cassist/cassist.h>
#endif
#include <pij/rank.h>
#ifdef	ENABLE_NETR
#include <netr/one.h>
#endif
#include <base/random.h>
#include <base/macros.h>
#include <base/lib.h>

#define MACROSTR(X)	#X
#define STR(X)	MACROSTR(X)

#define	SLIBINFONAME	STR(LIBINFONAME)
#define	SLIBINFOVERSION	STR(LIBINFOVERSION)
#define	LOGLV_AUTO	6

//Checks binary and library versions are the same. Otherwise prints warning
int checkversion()
{
	const char bln[]=SLIBINFONAME;
	const char blv[]=SLIBINFOVERSION;
	const char* ln=lib_name();
	const char* lv=lib_version();
	size_t lv1=lib_version1();
	size_t lv2=lib_version2();
	size_t lv3=lib_version3();
	int diff=0;
	
	if(strcmp(bln,ln))
	{
		diff=1;
		fprintf(stderr,"Error: different names for binary interface (%s) and library (%s). Skipped.%s",bln,ln,_NEWLINE_);
	}
	if(LIBINFOVERSION1!=lv1)
	{
		diff=1;
		fprintf(stderr,"Error: different major versions for binary interface (%s) and library (%s). Skipped.%s",blv,lv,_NEWLINE_);
	}
	else if(LIBINFOVERSION2!=lv2)
	{
		diff=1;
		fprintf(stderr,"Error: different minor versions for binary interface (%s) and library (%s). Skipped.%s",blv,lv,_NEWLINE_);
	}
	else if(LIBINFOVERSION3!=lv3)
		fprintf(stderr,"Warning: different minor versions for binary interface (%s) and library (%s). Loaded.%s",blv,lv,_NEWLINE_);
	if(!diff)
		fprintf(stderr,"Binary interface for library %s, version %s.%s",bln,blv,_NEWLINE_);
	return diff;
}

//Display usage info.
void usage(FILE* fp,const char* binpath)
{
	fprintf(fp,"NOTE: This is a brief usage instruction. For more details, see doc.pdf.%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth method [method-args...]%s",binpath,_NEWLINE_);
	fprintf(fp,"\tloglv:\tLog level, 0 for automatic (warnings + errors)%s",_NEWLINE_);
	fprintf(fp,"\trs:\tInitial random seed. 0 for using current time.%s",_NEWLINE_);
	fprintf(fp,"\tnth:\tMaximum number of parallel threads, 0 for automatic%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"For each method, look up the function of the same name in the library headers for usage.%s",_NEWLINE_);
	fprintf(fp,"Supports the following methods:%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp," Pairwise regulation probability (local precision) estimation%s",_NEWLINE_);
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"pij_gassist: invokes function pij_gassist in library and provides recommended causal inference of regulatory relations. \"gassist\" functions use discrete data as instrumental variable.%s",_NEWLINE_);
	fprintf(fp,"pij_gassist_trad: invokes function pij_gassist_trad in library and performs traditional causal inference test. WARNING: This is not and is not intended as a loyal reimplementation of Trigger R package.%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pij_gassist fg ft ft2 nt nt2 ns fp na nodiag memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Or:    %s loglv rs nth pij_gassist_trad fg ft ft2 nt nt2 ns fp na nodiag memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfg:\tPath of input file containing genotype data for E(A). Element [i,j] is the genotype value of the best eQTL of gene i of sample j, and should be among values 0,1,...,na. Data matrix size: (nt,ns). %s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp:\tPath of output file for estimated probability (local precision, e.g. 1-local FDR) of the test of interest, either the recommended test by pij_gassist, or the traditional test by pij_gassist_trad. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tna:\tNumber of allleles for the species considered.%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"pijs_gassist: invokes function pijs_gassist in library and performs five tests for causal inference test to estimate probabilities (local precisions) of pairwise regulation A->B with genotype information%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pijs_gassist fg ft ft2 nt nt2 ns fp1 fp2 fp3 fp4 fp5 na nodiag memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfg:\tPath of input file containing genotype data for E(A). Element [i,j] is the genotype value of the best eQTL of gene i of sample j, and should be among values 0,1,...,na. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp1:\tPath of output file for probability of test 1, E->A v.s. E  A. For nodiag=0, because the function expects significant eQTLs, p1 always return 1. For nodiag=1, uses diagonal elements of p2. Consider replacing p1 with your own (1-FDR) from eQTL discovery. Data vector size: (nt).%s",_NEWLINE_);
	fprintf(fp,"\tfp2:\tPath of output file for probability of test 2, E->B v.s. E  B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp3:\tPath of output file for probability of test 3, E->A->B v.s. E->A->B with E->B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp4:\tPath of output file for probability of test 4, E->A->B with E->B v.s. E->A  B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp5:\tPath of output file for probability of test 5, E->A->B with E->B v.s. A<-E->B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tna:\tNumber of allleles for the species considered.%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
#ifdef ENABLE_CASSIST
	fprintf(fp,"pij_cassist: invokes function pij_cassist in library and provides recommended causal inference of regulatory relations. \"cassist\" functions use continuous data as instrumental variable.%s",_NEWLINE_);
	fprintf(fp,"pij_cassist_trad: invokes function pij_cassist_trad in library and performs traditional causal inference test.%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pij_cassist fc ft ft2 nt nt2 ns fp nodiag memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Or:    %s loglv rs nth pij_cassist_trad fc ft ft2 nt nt2 ns fp nodiag memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfc:\tPath of input file containing continuous instrumental variable data for E(A). Element [i,j] is the value of the best continuous instrumental variable of gene i of sample j. Data matrix size: (nt,ns). %s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp:\tPath of output file for inferred posterior probability of the test of interest, either the recommended test by pij_gassist, or the traditional test by pij_gassist_trad. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"pijs_cassist: invokes function pijs_cassist in library and performs five tests for causal inference test to estimate probabilities (local precisions) of pairwise regulation A->B with continuous anchor information%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pijs_cassist fc ft ft2 nt nt2 ns fp1 fp2 fp3 fp4 fp5 nodiag memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfc:\tPath of input file containing continuous instrumental variable data for E(A). Element [i,j] is the value of the best continuous instrumental variable of gene i of sample j. Data matrix size: (nt,ns). %s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp1:\tPath of output file for probability of test 1, E->A v.s. E  A. For nodiag=0, because the function expects significant instrumental variables, p1 always return 1. For nodiag=1, uses diagonal elements of p2. Consider replacing p1 with your own (1-FDR) from instrumental variable discovery. Data vector size: (nt).%s",_NEWLINE_);
	fprintf(fp,"\tfp2:\tPath of output file for probability of test 2, E->B v.s. E  B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp3:\tPath of output file for probability of test 3, E->A->B v.s. E->A->B with E->B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp4:\tPath of output file for probability of test 4, E->A->B with E->B v.s. E->A  B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp5:\tPath of output file for probability of test 5, E->A->B with E->B v.s. A<-E->B. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
#endif
	fprintf(fp,"pij_rank: invokes function pij_rank in library to estimate probability (local precision) of pairwise correlation A--B with only pairwise expression data%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pij_rank ft ft2 nt nt2 ns fp nodiag memlimit%s",binpath, _NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp:\tPath of output file for probability. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tnodiag:\tWhether diagonal elements of output probability matrices should be neglected (due to identical A,B)%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage. %s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
#ifdef ENABLE_PV
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp," Pairwise regulation p-value estimation%s",_NEWLINE_);
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"pijs_gassist_pv: invokes function pijs_gassist_pv in library and performs five tests for causal inference test to estimate p-values of pairwise regulation A->B with genotype information%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pijs_gassist_pv fg ft ft2 nt nt2 ns fp1 fp2 fp3 fp4 fp5 na memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfg:\tPath of input file containing genotype data for E(A). Element [i,j] is the genotype value of the best eQTL of gene i of sample j, and should be among values 0,1,...,na. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp1:\tPath of output file for p-value of LLR of test 1, E->A (alt) v.s. E  A (null). Data vector size: (nt).%s",_NEWLINE_);
	fprintf(fp,"\tfp2:\tPath of output file for p-value of LLR of test 2, E->B (alt) v.s. E  B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp3:\tPath of output file for p-value of LLR of test 3, E->A->B with E->B (alt) v.s. E->A->B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp4:\tPath of output file for p-value of LLR of test 4, E->A->B with E->B (alt) v.s. E->A  B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp5:\tPath of output file for p-value of LLR of test 5, E->A->B with E->B (alt) v.s. A<-E->B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tna:\tNumber of allleles for the species considered.%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
#ifdef ENABLE_CASSIST
	fprintf(fp,"pijs_cassist_pv: invokes function pijs_cassist_pv in library and performs five tests for causal inference test to estimate p-values of pairwise regulation A->B with continuous anchor information%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pijs_cassist_pv fc ft ft2 nt nt2 ns fp1 fp2 fp3 fp4 fp5 memlimit%s",binpath,_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfc:\tPath of input file containing continuous instrumental variable data for E(A). Element [i,j] is the value of the best continuous instrumental variable of gene i of sample j. Data matrix size: (nt,ns). %s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp1:\tPath of output file for p-value of LLR of test 1, E->A (alt) v.s. E  A (null). Data vector size: (nt).%s",_NEWLINE_);
	fprintf(fp,"\tfp2:\tPath of output file for p-value of LLR of test 2, E->B (alt) v.s. E  B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp3:\tPath of output file for p-value of LLR of test 3, E->A->B with E->B (alt) v.s. E->A->B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp4:\tPath of output file for p-value of LLR of test 4, E->A->B with E->B (alt) v.s. E->A  B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tfp5:\tPath of output file for p-value of LLR of test 5, E->A->B with E->B (alt) v.s. A<-E->B (null). Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
#endif
	fprintf(fp,"pij_rank_pv: invokes function pij_rank_pv in library to estimate p-value of pairwise correlation A--B with only pairwise expression data%s",_NEWLINE_);
	fprintf(fp,"Null hypothesis: no correlation. Alternative hypothesis: allows correlation%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth pij_rank_pv ft ft2 nt nt2 ns fp memlimit%s",binpath, _NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tft:\tPath of input file containing gene expression data for A. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt,ns).%s",_NEWLINE_);
	fprintf(fp,"\tft2:\tPath of input file containing gene expression data for B. Element [i,j] is the expression level of gene i of sample j. Data matrix size: (nt2,ns).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of A's as integer%s",_NEWLINE_);
	fprintf(fp,"\tnt2:\tNumber of B's as integer%s",_NEWLINE_);
	fprintf(fp,"\tns:\tNumber of samples as integer%s",_NEWLINE_);
	fprintf(fp,"\tfp:\tPath of output file for p-value of LLR. Data matrix size: (nt,nt2).%s",_NEWLINE_);
	fprintf(fp,"\tmemlimit:\tThe approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage. %s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
#endif
#ifdef ENABLE_NETR
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp," Reconstruction of directed acyclic graph from prior information%s",_NEWLINE_);
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"netr_one_greedy: invokes function netr_one_greedy in library and reconstructs directed acyclic graph by adding one most significant edge at a time while avoiding cycles, according to the prior information on edge significance.%s",_NEWLINE_);
	fprintf(fp,"Usage: %s loglv rs nth netr_one_greedy fprior nt fnet namax nimax nomax%s",binpath,_NEWLINE_);
	fprintf(fp,"Method arguments:%s",_NEWLINE_);
	fprintf(fp,"\tfprior:\tPath of input file containing prior information of edge significance. Larger values indicate stronger significance. Example inputs are the pairwise regulation probability or minus pairwise regulation p-value outputs from above pij_* functions. Data matrix size: (nt,nt).%s",_NEWLINE_);
	fprintf(fp,"\tnt:\tNumber of nodes as integer%s",_NEWLINE_);
	fprintf(fp,"\tfnet:\tPath of output file for the reconstructed directed acyclic graph. Data matrix size: (nt,nt). [i,j]=1 indicates a directed edge from node i to node j, or [i,j]=0 otherwise. Data type: 8-bit unsigned integer.%s",_NEWLINE_);
	fprintf(fp,"\tnamax:\tThe maximum total number of edges for the reconstructed graph. Default 0 indicates unlimited (=nt*(nt-1)/2).%s",_NEWLINE_);
	fprintf(fp,"\tnimax:\tThe maximum number of incoming edges for each node. Default 0 indicates unlimited.%s",_NEWLINE_);
	fprintf(fp,"\tnomax:\tThe maximum number of outgoing edges for each node. Default 0 indicates unlimited.%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
#endif
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp," File formats:%s",_NEWLINE_);
	fprintf(fp,"**************************************************************%s",_NEWLINE_);
	fprintf(fp,"%s",_NEWLINE_);
	fprintf(fp,"This binary interface accepts two input/output file formats: raw and tsv. Within each run, all input and output files must all have the same format, either raw or tsv. This is indicated from the method name. Normal method names indicate raw file format, and appending the method name with '_tsv' would indicate tsv file format.%s",_NEWLINE_);
	fprintf(fp,"\tRaw format: for matrices, row-major sequence is used. All genotype data follow type 'GTYPE', and all expression data and output probabilities follow type 'FTYPE'. Detailed definitions of GTYPE and FTYPE can be found in base/types.h in source tree. By default, GTYPE is 8-bit (1-byte) unsigned char, FTYPE is 32-bit (4-byte) signed float in native endianness.%s",_NEWLINE_);
	fprintf(fp,"\tTsv format: tab separated value files. Each row should be separated by new line, and each column by tab or space. The number of rows and columns must match input dimensional parameters.%s",_NEWLINE_);
}

//vector and matrix input output functions
static int bin_checkfilesize(FILE* fp,size_t s,const char f[])
{
	long fsize;
	if(fseek(fp,0,SEEK_END))
	{
		LOG(1,"Can't fseek file %s.",f)
		return 1;
	}
	fsize=ftell(fp);
	if(fsize<0)
	{
		LOG(1,"Can't ftell file %s.",f)
		return 1;
	}
	if((size_t)fsize!=s)
	{
		LOG(1,"Incorrect file size for %s",f)
		return 1;
	}
	if(fseek(fp,0,SEEK_SET))
	{
		LOG(1,"Can't fseek file %s.",f)
		return 1;
	}
	return 0;
}

#define	OFORMAT	"%.8G"

static int bin_matrixfi_raw(const char f[],MATRIXF* m)
{
	FILE	*fp;
	fp=fopen(f,"rb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(bin_checkfilesize(fp,m->size1*m->size2*sizeof(m->data[0]),f))
	{
		fclose(fp);
		return 1;
	}
	if(MATRIXFF(fread)(fp,m))
	{
		LOG(1,"Can't read file or has wrong size: %s.",f)
		fclose(fp);
		return 1;
	}
	fclose(fp);	
	return 0;
}

static int bin_matrixfi_tsv(const char f[],MATRIXF* m)
{
	FILE	*fp;
	FTYPE	tf;
	fp=fopen(f,"rb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(MATRIXFF(fscanf)(fp,m))
	{
		LOG(1,"Can't read file or has wrong size: %s.",f)
		fclose(fp);
		return 1;
	}
	if(fscanf(fp,"%f",&tf)==1)
	{
		LOG(1,"File has more data entries than expected: %s.",f)
		fclose(fp);
		return 1;
	}
	fclose(fp);	
	return 0;
}

static int bin_matrixgi_raw(const char f[],MATRIXG* m)
{
	FILE	*fp;
	fp=fopen(f,"rb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(bin_checkfilesize(fp,m->size1*m->size2*sizeof(m->data[0]),f))
	{
		fclose(fp);
		return 1;
	}
	if(MATRIXGF(fread)(fp,m))
	{
		LOG(1,"Can't read file or has wrong size: %s.",f)
		fclose(fp);
		return 1;
	}

	fclose(fp);	
	return 0;
}

static int bin_matrixgi_tsv(const char f[],MATRIXG* m)
{
	FILE	*fp;
	int		tg;
	fp=fopen(f,"rb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(MATRIXGF(fscanf)(fp,m))
	{
		LOG(1,"Can't read file or has wrong size: %s.",f)
		fclose(fp);
		return 1;
	}
	if(fscanf(fp,"%i",&tg)==1)
	{
		LOG(1,"File has more data entries than expected: %s.",f)
		fclose(fp);
		return 1;
	}
	fclose(fp);	
	return 0;
}

static int bin_matrixfo_raw(const char f[],const MATRIXF* m)
{
	FILE	*fp;
	fp=fopen(f,"wb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(MATRIXFF(fwrite)(fp,m))
	{
		LOG(1,"Can't write file: %s.",f)
		fclose(fp);
		return 1;
	}
	fclose(fp);	
	return 0;
}

static int bin_matrixfo_tsv(const char f[],const MATRIXF* m)
{
	FILE	*fp;
	size_t	i,j;
	fp=fopen(f,"wb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	for(j=0;j<m->size1;j++)
	{
		if(!fprintf(fp,OFORMAT,MATRIXFF(get)(m,j,0)))
		{
			LOG(1,"Can't write file: %s.",f)
			fclose(fp);
			return 1;
		}
		for(i=1;i<m->size2;i++)
			if((!putc('\t',fp))||(!fprintf(fp,OFORMAT,MATRIXFF(get)(m,j,i))))
			{
				LOG(1,"Can't write file: %s.",f)
				fclose(fp);
				return 1;
			}
		if(!fputs(_NEWLINE_,fp))
		{
			LOG(1,"Can't write file: %s.",f)
			fclose(fp);
			return 1;
		}
	}
	fclose(fp);	
	return 0;
}

static int bin_matrixuco_raw(const char f[],const MATRIXUC* m)
{
	FILE	*fp;
	fp=fopen(f,"wb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(MATRIXUCF(fwrite)(fp,m))
	{
		LOG(1,"Can't write file: %s.",f)
		fclose(fp);
		return 1;
	}
	fclose(fp);	
	return 0;
}

static int bin_matrixuco_tsv(const char f[],const MATRIXUC* m)
{
	FILE	*fp;
	size_t	i,j;
	fp=fopen(f,"wb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	for(j=0;j<m->size1;j++)
	{
		if(!fprintf(fp,"%u",MATRIXUCF(get)(m,j,0)))
		{
			LOG(1,"Can't write file: %s.",f)
			fclose(fp);
			return 1;
		}
		for(i=1;i<m->size2;i++)
			if((!putc('\t',fp))||(!fprintf(fp,"%u",MATRIXUCF(get)(m,j,i))))
			{
				LOG(1,"Can't write file: %s.",f)
				fclose(fp);
				return 1;
			}
		if(!fputs(_NEWLINE_,fp))
		{
			LOG(1,"Can't write file: %s.",f)
			fclose(fp);
			return 1;
		}
	}
	fclose(fp);	
	return 0;
}

static int bin_vectorfo_raw(const char f[],const VECTORF* v)
{
	FILE	*fp;
	fp=fopen(f,"wb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(VECTORFF(fwrite)(fp,v))
	{
		LOG(1,"Can't write file: %s.",f)
		fclose(fp);
		return 1;
	}
	fclose(fp);	
	return 0;
}

static int bin_vectorfo_tsv(const char f[],const VECTORF* v)
{
	FILE	*fp;
	size_t	i;
	fp=fopen(f,"wb");
	if(!fp)
	{
		LOG(1,"Can't open file %s.",f)
		return 1;
	}
	if(!fprintf(fp,OFORMAT,VECTORFF(get)(v,0)))
	{
		LOG(1,"Can't write file: %s.",f)
		fclose(fp);
		return 1;
	}
	for(i=1;i<v->size;i++)
		if((!putc('\t',fp))||(!fprintf(fp,OFORMAT,VECTORFF(get)(v,i))))
		{
			LOG(1,"Can't write file: %s.",f)
			fclose(fp);
			return 1;
		}
	if(!fputs(_NEWLINE_,fp))
	{
		LOG(1,"Can't write file: %s.",f)
		fclose(fp);
		return 1;
	}
	fclose(fp);
	return 0;
}
#undef	OFORMAT

//Library initialization
int bin_call_lib_init(const char* argv[])
{
	unsigned char loglv;
	unsigned long rs;
	size_t nth;
	int	loglv0,nth0;
	loglv0=(unsigned char)atoi(argv[0]);
	if((loglv0<0)||(loglv0>12))
	{
		fprintf(stderr,"Log level must be between 0 and 12.");
		return 1;
	}
	loglv=(unsigned char)loglv0;
	if(!loglv)
		loglv=LOGLV_AUTO;
	sscanf(argv[1],"%lu",&rs);
	nth0=atoi(argv[2]);
	if(nth0<0)
	{
		fprintf(stderr,"Thead count must be nonnegative.");
		return 1;
	}
	nth=(size_t)nth0;

	lib_init(loglv,rs,nth);
	if(nth>64)
		LOG(6,"Thread count %lu is very large.", nth)
	return 0;
}

/*****************************************************************************
 *	Exposed functions
 ****************************************************************************/

//pij_rank generic function, for raw and tsv versions
int bin_pij_rank_func(int argc,const char* argv[],int (*fin_fm)(const char[],MATRIXF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATF(t)CLEANMATF(t2)CLEANMATF(p)
	const char *f_t,*f_t2,*f_p;
	char nodiag;
	size_t	ns,ng,nt,memlimit;
	const char *funcname=argv[0];
	MATRIXF	*t=0,*t2=0,*p=0;
	
	if(argc!=9)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_t=argv[1];
	f_t2=argv[2];
	ng=(size_t)atol(argv[3]);
	nt=(size_t)atol(argv[4]);
	ns=(size_t)atol(argv[5]);
	if(!(ng&&nt&&ns))
	{
		LOG(0,"Invalid input dimensions.")
		return -1;
	}
	f_p=argv[6];
	nodiag=(char)atoi(argv[7]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value %i.",nodiag)
		return -1;
	}
	memlimit=(size_t)atol(argv[8]);
	if(!memlimit)
		memlimit=(size_t)-1;
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	t=MATRIXFF(alloc)(ng,ns);
	t2=MATRIXFF(alloc)(nt,ns);
	p=MATRIXFF(alloc)(ng,nt);
	if(!(t&&t2&&p))
		ERRRET("Not enough memory.")
	
	//File reads
	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(pij_rank(t,t2,p,nodiag,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	LOG(11,"Writing file %s.",f_p)
	if(fout_fm(f_p,p))
		ERRRET("Can't write to file %s.",f_p)

	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_pij_rank(int argc,const char* argv[])
{
	return bin_pij_rank_func(argc,argv,bin_matrixfi_raw,bin_matrixfo_raw);
}

static inline int bin_pij_rank_tsv(int argc,const char* argv[])
{
	return bin_pij_rank_func(argc,argv,bin_matrixfi_tsv,bin_matrixfo_tsv);
}

//pij_gassist generic function, for pij_gassist and pij_gassist_trad, for raw and tsv versions
int bin_pij_gassist_func(int argc,const char* argv[],int (*func)(const MATRIXG*,const MATRIXF*,const MATRIXF*,MATRIXF*,size_t,char,size_t),int (*fin_gm)(const char[],MATRIXG*),int (*fin_fm)(const char[],MATRIXF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATG(g)CLEANMATF(t)CLEANMATF(t2)CLEANMATF(p)
	const char *f_g,*f_t,*f_t2,*f_p;
	char nodiag;
	const char *funcname=argv[0];
	size_t	ns,ng,nt,memlimit;
	size_t	nv;
	MATRIXG		*g=0;
	MATRIXF		*t=0,*t2=0,*p=0;
	
	if(argc!=11)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[1];
	f_t=argv[2];
	f_t2=argv[3];
	ng=(size_t)atol(argv[4]);
	nt=(size_t)atol(argv[5]);
	ns=(size_t)atol(argv[6]);
	f_p=argv[7];
	nv=(size_t)atol(argv[8])+1;
	if((!(ng&&nt&&ns))||!nv)
	{
		LOG(0,"Invalid input dimensions or allele count.")
		return -1;
	}
	nodiag=(char)atoi(argv[9]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value.")
		return -1;
	}
	memlimit=(size_t)atol(argv[10]);
	if(!memlimit)
		memlimit=(size_t)-1;
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
	if(fin_gm(f_g,g))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_g)
	if(nv==1)
	{
		nv=(size_t)MATRIXGF(max)(g)+1;
		if(nv<2)
			ERRRET("Failed to autodetect allele count. Please check your genotype input file.")
		else
			LOG(9,"Autodetected number of alleles=%u.",nv-1)
	}

	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(func(g,t,t2,p,nv,nodiag,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	if(fout_fm(f_p,p))
		ERRRET("Can't write file %s.",f_p)
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef FWRITE
#undef CLEANUP
}

static inline int bin_pij_gassist(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist,bin_matrixgi_raw,bin_matrixfi_raw,bin_matrixfo_raw);
}

static inline int bin_pij_gassist_tsv(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist,bin_matrixgi_tsv,bin_matrixfi_tsv,bin_matrixfo_tsv);
}

static inline int bin_pij_gassist_trad(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist_trad,bin_matrixgi_raw,bin_matrixfi_raw,bin_matrixfo_raw);
}

static inline int bin_pij_gassist_trad_tsv(int argc,const char* argv[])
{
	return bin_pij_gassist_func(argc,argv,pij_gassist_trad,bin_matrixgi_tsv,bin_matrixfi_tsv,bin_matrixfo_tsv);
}

//pijs_gassist generic function, for pijs_gassist, for raw and tsv versions
int bin_pijs_gassist_func(int argc,const char* argv[],int (*fin_gm)(const char[],MATRIXG*),int (*fin_fm)(const char[],MATRIXF*),int (*fout_fv)(const char[],const VECTORF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATG(g)CLEANMATF(t)CLEANMATF(t2)CLEANVECF(p1)\
				CLEANMATF(p2)CLEANMATF(p3)CLEANMATF(p4)CLEANMATF(p5)
	const char *f_g,*f_t,*f_t2,*f_p1,*f_p2,*f_p3,*f_p4,*f_p5;
	char nodiag;
	size_t	ns,ng,nt,memlimit;
	size_t	nv;
	const char *funcname=argv[0];
	VECTORF		*p1=0;
	MATRIXG		*g=0;
	MATRIXF		*t=0,*t2=0,*p2=0,*p3=0,*p4=0,*p5=0;
	
	if(argc!=15)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[1];
	f_t=argv[2];
	f_t2=argv[3];
	ng=(size_t)atol(argv[4]);
	nt=(size_t)atol(argv[5]);
	ns=(size_t)atol(argv[6]);
	f_p1=argv[7];
	f_p2=argv[8];
	f_p3=argv[9];
	f_p4=argv[10];
	f_p5=argv[11];
	nv=(size_t)atol(argv[12])+1;
	if((!(ng&&nt&&ns))||!nv)
	{
		LOG(0,"Invalid input dimensions or allele count.")
		return -1;
	}
	nodiag=(char)atoi(argv[13]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value.")
		return -1;
	}
	memlimit=(size_t)atol(argv[14]);
	if(!memlimit)
		memlimit=(size_t)-1;
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
	if(fin_gm(f_g,g))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_g)
	if(nv==1)
	{
		nv=(size_t)MATRIXGF(max)(g)+1;
		if(nv<2)
			ERRRET("Failed to autodetect allele count. Please check your genotype input file.")
		else
			LOG(9,"Autodetected number of alleles=%u.",nv-1)
	}
	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(pijs_gassist(g,t,t2,p1,p2,p3,p4,p5,nv,nodiag,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	if(fout_fv(f_p1,p1))
		ERRRET("Can't write file %s.",f_p1)
	if(fout_fm(f_p2,p2))
		ERRRET("Can't write file %s.",f_p2)
	if(fout_fm(f_p3,p3))
		ERRRET("Can't write file %s.",f_p3)
	if(fout_fm(f_p4,p4))
		ERRRET("Can't write file %s.",f_p4)
	if(fout_fm(f_p5,p5))
		ERRRET("Can't write file %s.",f_p5)
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_pijs_gassist(int argc,const char* argv[])
{
	return bin_pijs_gassist_func(argc,argv,bin_matrixgi_raw,bin_matrixfi_raw,bin_vectorfo_raw,bin_matrixfo_raw);
}

static inline int bin_pijs_gassist_tsv(int argc,const char* argv[])
{
	return bin_pijs_gassist_func(argc,argv,bin_matrixgi_tsv,bin_matrixfi_tsv,bin_vectorfo_tsv,bin_matrixfo_tsv);
}

#ifdef	ENABLE_CASSIST
//pijs_cassist generic function, for pijs_cassist, for raw and tsv versions
int bin_pijs_cassist_func(int argc,const char* argv[],int (*fin_fm)(const char[],MATRIXF*),int (*fout_fv)(const char[],const VECTORF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATF(g)CLEANMATF(t)CLEANMATF(t2)CLEANVECF(p1)\
				CLEANMATF(p2)CLEANMATF(p3)CLEANMATF(p4)CLEANMATF(p5)
	const char *f_g,*f_t,*f_t2,*f_p1,*f_p2,*f_p3,*f_p4,*f_p5;
	char nodiag;
	const char *funcname=argv[0];
	size_t	ns,ng,nt,memlimit;
	VECTORF		*p1=0;
	MATRIXF		*g=0;
	MATRIXF		*t=0,*t2=0,*p2=0,*p3=0,*p4=0,*p5=0;
	
	if(argc!=14)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[1];
	f_t=argv[2];
	f_t2=argv[3];
	ng=(size_t)atol(argv[4]);
	nt=(size_t)atol(argv[5]);
	ns=(size_t)atol(argv[6]);
	f_p1=argv[7];
	f_p2=argv[8];
	f_p3=argv[9];
	f_p4=argv[10];
	f_p5=argv[11];
	if(!(ng&&nt&&ns))
	{
		LOG(0,"Invalid input dimensions.")
		return -1;
	}
	nodiag=(char)atoi(argv[12]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value.")
		return -1;
	}
	memlimit=(size_t)atol(argv[13]);
	if(!memlimit)
		memlimit=(size_t)-1;
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	g=MATRIXFF(alloc)(ng,ns);
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
	if(fin_fm(f_g,g))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_g)
	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(pijs_cassist(g,t,t2,p1,p2,p3,p4,p5,nodiag,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	if(fout_fv(f_p1,p1))
		ERRRET("Can't write file %s.",f_p1)
	if(fout_fm(f_p2,p2))
		ERRRET("Can't write file %s.",f_p2)
	if(fout_fm(f_p3,p3))
		ERRRET("Can't write file %s.",f_p3)
	if(fout_fm(f_p4,p4))
		ERRRET("Can't write file %s.",f_p4)
	if(fout_fm(f_p5,p5))
		ERRRET("Can't write file %s.",f_p5)
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_pijs_cassist(int argc,const char* argv[])
{
	return bin_pijs_cassist_func(argc,argv,bin_matrixfi_raw,bin_vectorfo_raw,bin_matrixfo_raw);
}

static inline int bin_pijs_cassist_tsv(int argc,const char* argv[])
{
	return bin_pijs_cassist_func(argc,argv,bin_matrixfi_tsv,bin_vectorfo_tsv,bin_matrixfo_tsv);
}

//pij_cassist generic function, for pij_cassist and pij_cassist_trad, for raw and tsv versions
int bin_pij_cassist_func(int argc,const char* argv[],int (*func)(const MATRIXF*,const MATRIXF*,const MATRIXF*,MATRIXF*,char,size_t),int (*fin_fm)(const char[],MATRIXF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATF(g)CLEANMATF(t)CLEANMATF(t2)CLEANMATF(p)
	const char *f_g,*f_t,*f_t2,*f_p;
	const char *funcname=argv[0];
	char nodiag;
	size_t	ns,ng,nt,memlimit;
	MATRIXF		*g=0;
	MATRIXF		*t=0,*t2=0,*p=0;
	
	if(argc!=10)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[1];
	f_t=argv[2];
	f_t2=argv[3];
	ng=(size_t)atol(argv[4]);
	nt=(size_t)atol(argv[5]);
	ns=(size_t)atol(argv[6]);
	f_p=argv[7];
	if(!(ng&&nt&&ns))
	{
		LOG(0,"Invalid input dimensions.")
		return -1;
	}
	nodiag=(char)atoi(argv[8]);
	if((nodiag!=0)&&(nodiag!=1))
	{
		LOG(0,"Invalid nodiag value.")
		return -1;
	}
	memlimit=(size_t)atol(argv[9]);
	if(!memlimit)
		memlimit=(size_t)-1;
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	g=MATRIXFF(alloc)(ng,ns);
	t=MATRIXFF(alloc)(ng,ns);
	t2=MATRIXFF(alloc)(nt,ns);
	p=MATRIXFF(alloc)(ng,nt);
	if(!(g&&t&&t2&&p))
		ERRRET("Not enough memory.")
	
	//File reads
	LOG(11,"Reading file %s.",f_g)
	if(fin_fm(f_g,g))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_g)
	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(func(g,t,t2,p,nodiag,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	if(fout_fm(f_p,p))
		ERRRET("Can't write file %s.",f_p)
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef FWRITE
#undef CLEANUP
}

static inline int bin_pij_cassist(int argc,const char* argv[])
{
	return bin_pij_cassist_func(argc,argv,pij_cassist,bin_matrixfi_raw,bin_matrixfo_raw);
}

static inline int bin_pij_cassist_tsv(int argc,const char* argv[])
{
	return bin_pij_cassist_func(argc,argv,pij_cassist,bin_matrixfi_tsv,bin_matrixfo_tsv);
}

static inline int bin_pij_cassist_trad(int argc,const char* argv[])
{
	return bin_pij_cassist_func(argc,argv,pij_cassist_trad,bin_matrixfi_raw,bin_matrixfo_raw);
}

static inline int bin_pij_cassist_trad_tsv(int argc,const char* argv[])
{
	return bin_pij_cassist_func(argc,argv,pij_cassist_trad,bin_matrixfi_tsv,bin_matrixfo_tsv);
}
#endif

#ifdef ENABLE_PV
//pij_rank_pv generic function, for raw and tsv versions
static int bin_pij_rank_pv_func(int argc,const char* argv[],int (*fin_fm)(const char[],MATRIXF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATF(t)CLEANMATF(t2)CLEANMATF(p)
	const char *f_t,*f_t2,*f_p;
	const char *funcname=argv[0];

	size_t	ns,ng,nt,memlimit;
	MATRIXF	*t=0,*t2=0,*p=0;
	
	if(argc!=8)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_t=argv[1];
	f_t2=argv[2];
	ng=(size_t)atol(argv[3]);
	nt=(size_t)atol(argv[4]);
	ns=(size_t)atol(argv[5]);
	if(!(ng&&nt&&ns))
	{
		LOG(0,"Invalid input dimensions.")
		return -1;
	}
	f_p=argv[6];
	memlimit=(size_t)atol(argv[7]);
	if(!memlimit)
		memlimit=(size_t)-1;
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	t=MATRIXFF(alloc)(ng,ns);
	t2=MATRIXFF(alloc)(nt,ns);
	p=MATRIXFF(alloc)(ng,nt);
	if(!(t&&t2&&p))
		ERRRET("Not enough memory.")
	
	//File reads
	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(pij_rank_pv(t,t2,p,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	LOG(11,"Writing file %s.",f_p)
	if(fout_fm(f_p,p))
		ERRRET("Can't write to file %s.",f_p)

	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_pij_rank_pv(int argc,const char* argv[])
{
	return bin_pij_rank_pv_func(argc,argv,bin_matrixfi_raw,bin_matrixfo_raw);
}

static inline int bin_pij_rank_pv_tsv(int argc,const char* argv[])
{
	return bin_pij_rank_pv_func(argc,argv,bin_matrixfi_tsv,bin_matrixfo_tsv);
}

//pijs_gassist_pv generic function, for pijs_gassist_pv, for raw and tsv versions
int bin_pijs_gassist_pv_func(int argc,const char* argv[],int (*fin_gm)(const char[],MATRIXG*),int (*fin_fm)(const char[],MATRIXF*),int (*fout_fv)(const char[],const VECTORF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATG(g)CLEANMATF(t)CLEANMATF(t2)CLEANVECF(p1)\
				CLEANMATF(p2)CLEANMATF(p3)CLEANMATF(p4)CLEANMATF(p5)
	const char *f_g,*f_t,*f_t2,*f_p1,*f_p2,*f_p3,*f_p4,*f_p5;
	const char *funcname=argv[0];
	size_t	ns,ng,nt,memlimit;
	size_t	nv;
	VECTORF		*p1=0;
	MATRIXG		*g=0;
	MATRIXF		*t=0,*t2=0,*p2=0,*p3=0,*p4=0,*p5=0;
	
	if(argc!=14)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[1];
	f_t=argv[2];
	f_t2=argv[3];
	ng=(size_t)atol(argv[4]);
	nt=(size_t)atol(argv[5]);
	ns=(size_t)atol(argv[6]);
	f_p1=argv[7];
	f_p2=argv[8];
	f_p3=argv[9];
	f_p4=argv[10];
	f_p5=argv[11];
	nv=(size_t)atol(argv[12])+1;
	if((!(ng&&nt&&ns))||!nv)
	{
		LOG(0,"Invalid input dimensions or allele count.")
		return -1;
	}
	memlimit=(size_t)atol(argv[13]);
	if(!memlimit)
		memlimit=(size_t)-1;
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
	if(fin_gm(f_g,g))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_g)
	if(nv==1)
	{
		nv=(size_t)MATRIXGF(max)(g)+1;
		if(nv<2)
			ERRRET("Failed to autodetect allele count. Please check your genotype input file.")
		else
			LOG(9,"Autodetected number of alleles=%u.",nv-1)
	}
	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(pijs_gassist_pv(g,t,t2,p1,p2,p3,p4,p5,nv,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	if(fout_fv(f_p1,p1))
		ERRRET("Can't write file %s.",f_p1)
	if(fout_fm(f_p2,p2))
		ERRRET("Can't write file %s.",f_p2)
	if(fout_fm(f_p3,p3))
		ERRRET("Can't write file %s.",f_p3)
	if(fout_fm(f_p4,p4))
		ERRRET("Can't write file %s.",f_p4)
	if(fout_fm(f_p5,p5))
		ERRRET("Can't write file %s.",f_p5)
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_pijs_gassist_pv(int argc,const char* argv[])
{
	return bin_pijs_gassist_pv_func(argc,argv,bin_matrixgi_raw,bin_matrixfi_raw,bin_vectorfo_raw,bin_matrixfo_raw);
}

static inline int bin_pijs_gassist_pv_tsv(int argc,const char* argv[])
{
	return bin_pijs_gassist_pv_func(argc,argv,bin_matrixgi_tsv,bin_matrixfi_tsv,bin_vectorfo_tsv,bin_matrixfo_tsv);
}

#ifdef ENABLE_CASSIST
//pijs_cassist_pv generic function, for pijs_cassist_pv, for raw and tsv versions
int bin_pijs_cassist_pv_func(int argc,const char* argv[],int (*fin_fm)(const char[],MATRIXF*),int (*fout_fv)(const char[],const VECTORF*),int (*fout_fm)(const char[],const MATRIXF*))
{
#define CLEANUP CLEANMATF(g)CLEANMATF(t)CLEANMATF(t2)CLEANVECF(p1)\
				CLEANMATF(p2)CLEANMATF(p3)CLEANMATF(p4)CLEANMATF(p5)
	const char *f_g,*f_t,*f_t2,*f_p1,*f_p2,*f_p3,*f_p4,*f_p5;
	const char *funcname=argv[0];
	size_t	ns,ng,nt,memlimit;
	VECTORF		*p1=0;
	MATRIXF		*g=0;
	MATRIXF		*t=0,*t2=0,*p2=0,*p3=0,*p4=0,*p5=0;
	
	if(argc!=13)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_g=argv[1];
	f_t=argv[2];
	f_t2=argv[3];
	ng=(size_t)atol(argv[4]);
	nt=(size_t)atol(argv[5]);
	ns=(size_t)atol(argv[6]);
	f_p1=argv[7];
	f_p2=argv[8];
	f_p3=argv[9];
	f_p4=argv[10];
	f_p5=argv[11];
	if((!(ng&&nt&&ns)))
	{
		LOG(0,"Invalid input dimensions.")
		return -1;
	}
	memlimit=(size_t)atol(argv[12]);
	if(!memlimit)
		memlimit=(size_t)-1;
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	g=MATRIXFF(alloc)(ng,ns);
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
	if(fin_fm(f_g,g))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_g)
	LOG(11,"Reading file %s.",f_t)
	if(fin_fm(f_t,t))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t)
	LOG(11,"Reading file %s.",f_t2)
	if(fin_fm(f_t2,t2))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_t2)
	
	//Calculation
	if(pijs_cassist_pv(g,t,t2,p1,p2,p3,p4,p5,memlimit))
		ERRRET("%s failed.",funcname)
	
	//File writes
	if(fout_fv(f_p1,p1))
		ERRRET("Can't write file %s.",f_p1)
	if(fout_fm(f_p2,p2))
		ERRRET("Can't write file %s.",f_p2)
	if(fout_fm(f_p3,p3))
		ERRRET("Can't write file %s.",f_p3)
	if(fout_fm(f_p4,p4))
		ERRRET("Can't write file %s.",f_p4)
	if(fout_fm(f_p5,p5))
		ERRRET("Can't write file %s.",f_p5)
	
	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_pijs_cassist_pv(int argc,const char* argv[])
{
	return bin_pijs_cassist_pv_func(argc,argv,bin_matrixfi_raw,bin_vectorfo_raw,bin_matrixfo_raw);
}

static inline int bin_pijs_cassist_pv_tsv(int argc,const char* argv[])
{
	return bin_pijs_cassist_pv_func(argc,argv,bin_matrixfi_tsv,bin_vectorfo_tsv,bin_matrixfo_tsv);
}

#endif
#endif

#ifdef ENABLE_NETR
//netr_one_greedy generic function, for raw and tsv versions
int bin_netr_one_greedy_func(int argc,const char* argv[],int (*fin_fm)(const char[],MATRIXF*),int (*fout_fm)(const char[],const MATRIXUC*))
{
#define CLEANUP CLEANMATF(p)CLEANMATUC(net)
	const char *f_i,*f_o;
	size_t	nt,namax,nimax,nomax,ret;
	long namax0,nimax0,nomax0;
	const char *funcname=argv[0];
	MATRIXF	*p=0;
	MATRIXUC	*net=0;
	
	if(argc!=7)
	{
		LOG(0,"Wrong argument count.")
		return -1;
	}
	f_i=argv[1];
	nt=(size_t)atol(argv[2]);
	if(!nt)
	{
		LOG(0,"Invalid input dimension.")
		return -1;
	}
	f_o=argv[3];
	namax0=atol(argv[4]);
	nimax0=atol(argv[5]);
	nomax0=atol(argv[6]);
	if((namax0<0)||(nimax0<0)||(nomax0<0))
	{
		LOG(0,"Edge count constraints must be non-negative. Use zero for unlimited.")
		return -1;
	}
	namax=(size_t)namax0;
	nimax=(size_t)nimax0;
	nomax=(size_t)nomax0;
	if(!namax)
		namax=(size_t)-1;
	if(!nimax)
		nimax=(size_t)-1;
	if(!nomax)
		nomax=(size_t)-1;
	LOG(8,"%s started.",funcname)
	
	//Memory allocation
	p=MATRIXFF(alloc)(nt,nt);
	net=MATRIXUCF(alloc)(nt,nt);
	if(!(p&&net))
		ERRRET("Not enough memory.")
	
	//File reads
	LOG(11,"Reading file %s.",f_i)
	if(fin_fm(f_i,p))
		ERRRET("Can't read file or has wrong size: %s. Make sure your file format matches with method name. For text or tsv files, use the _tsv suffix.",f_i)

	//Calculation
	if(!(ret=netr_one_greedy(p,net,namax,nimax,nomax)))
		ERRRET("Construction of a single network with netr_one_greedy failed.")
	else
		LOG(10,"Network constructed with a total of %lu edges.",ret)
	
	//File writes
	LOG(11,"Writing file %s.",f_o)
	if(fout_fm(f_o,net))
		ERRRET("Can't write to file %s.",f_o)

	CLEANUP
	LOG(8,"%s completed.",funcname)
	return 0;
#undef CLEANUP
}

static inline int bin_netr_one_greedy(int argc,const char* argv[])
{
	return bin_netr_one_greedy_func(argc,argv,bin_matrixfi_raw,bin_matrixuco_raw);
}

static inline int bin_netr_one_greedy_tsv(int argc,const char* argv[])
{
	return bin_netr_one_greedy_func(argc,argv,bin_matrixfi_tsv,bin_matrixuco_tsv);
}

#endif

int main(int argc,const char* argv[])
{
#define	GETFUNCNAMERAW(X)	if(!strcmp(argv[4],STR(X))) func=bin_##X;
#define	GETFUNCNAMETSV(X)	if(!strcmp(argv[4],STR(X)"_tsv")) func=bin_##X##_tsv;
#define	GETFUNCNAMECSV(X)	if(!strcmp(argv[4],STR(X)"_csv")) func=bin_##X##_tsv;
#define	GETFUNCNAMES(X)		GETFUNCNAMERAW(X) else GETFUNCNAMETSV(X) \
							else GETFUNCNAMECSV(X)
	int	ret;
	int	(*func)(int,const char*[]);
	
	if(checkversion())
		return 1;
	func=0;
	//Look for method name
	if(argc>=5)
	{
		GETFUNCNAMES(pijs_gassist)
		else GETFUNCNAMES(pij_gassist)
		else GETFUNCNAMES(pij_gassist_trad)
#ifdef ENABLE_CASSIST
		else GETFUNCNAMES(pijs_cassist)
		else GETFUNCNAMES(pij_cassist)
		else GETFUNCNAMES(pij_cassist_trad)
#endif
		else GETFUNCNAMES(pij_rank)
#ifdef ENABLE_PV
		else GETFUNCNAMES(pijs_gassist_pv)
	#ifdef ENABLE_CASSIST
		else GETFUNCNAMES(pijs_cassist_pv)
	#endif
		else GETFUNCNAMES(pij_rank_pv)
#endif
#ifdef ENABLE_NETR
		else GETFUNCNAMES(netr_one_greedy)
#endif
	}
	//If method name not found
	if(!func)
	{
		usage(stdout,argv[0]);
		return (argc>1);
	}
	if((ret=bin_call_lib_init(argv+1)))
		return ret;
	ret=func(argc-4,argv+4);
	return ret;
}
#undef	GETFUNCNAME
#undef	LIBINFO


































