// Main header file
// Author: Yi Wang (yi dot wang at unsw dot edu dot au
// 16-Jun-2011

#ifndef _RESAMPTEST_H
#define _RESAMPTEST_H

//#define MATHLIB_STANDALONE 
//#include "/usr/local/R/2.13/lib64/R/include/Rmath.h"
#include "Rmath.h"
#include "R.h"
/*//Rmath
#define Rf_dpois dpois
#define Rf_ppois ppois
#define Rf_qpois qpois
#define Rf_rpois rpois
#define Rf_dbinom dbinom
#define Rf_pbinom pbinom
#define Rf_qbinom qbinom 
#define Rf_rbinom rbinom 
#define Rf_dnbinom dnbinom
#define Rf_pnbinom pnbinom
#define Rf_qnbinom qnbinom
#define Rf_rnbinom rnbinom
*/

#define printf Rprintf

#include <inttypes.h> // To use the uintptr_t and intptr_t to make code portable
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_errno.h>

// rmv.h
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_gamma.h>

// return status
#define SUCCESS 0
#define FAILED 1
#define CannotOpenFile 2  
// logic
#define TRUE 1
#define FALSE 0
// model
#define LM 0
#define POISSON 1
#define NB 2
#define BIN 3
// shrinkage
#define NOSHRINK 0
#define IDENTITY 1
#define SHRINK 2
// test
#define LOGWILK 0
#define HOTELING 1
#define WALD 2
#define SCORE 3
#define LR 4
// estiMethod
#define NEWTON 0
#define CHI2 1
#define FISHER 2
// infoMatrix
#define OIM 0
#define EIM 1
// resampling
#define CASEBOOT 0
#define RESIBOOT 1
#define SCOREBOOT 2
#define PERMUTE 3
#define FREEPERM 4
#define MONTECARLO 5
#define RESIZ 6 
#define SCOREZ 7
#define PITSBOOT 8
// p-value adjustment
#define NONE 0
#define UNADJUST 1
#define FREESTEP 2
#define SINGLESTEP 3
#define STEPUP 4
#define NOMONO 5
// R-squared
#define HOOPER 0
#define VECTOR 1
// others
#define TOL 1e-6
#define MAXITER 999 
#define LAMBDA 0.8  // no shrinkage 
#define NaN -100000
#define MAX_LINE_LENGTH 65536
#define WRAP 4
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))
#define ABS(a) GSL_MAX(a, -a)  
//simulations
#define NSIMU 500 
#define NQ 6
#define RHO 0.5
#define ALFA 0.05
#define MAXLINELEN 256

typedef struct MethodStruc {
    // hypo test methods
    uintptr_t nboot;
    uintptr_t corr;
    uintptr_t test;
    uintptr_t resamp;
    uintptr_t reprand;
    uintptr_t student;
    uintptr_t punit;
    uintptr_t rsquare;
    uintptr_t nRows;
    uintptr_t nVars;
    uintptr_t nParam;
    uintptr_t showtime;

    // numeric
    double shrink_param;
    gsl_vector *smry_lambda;
    gsl_vector *anova_lambda;
    double tol;
} mv_Method;

// used for manylm only
typedef struct matStruc {
    gsl_matrix *mat;   // hat(X)
    gsl_matrix *SS;
    gsl_matrix *Coef;
    gsl_matrix *Res;
    gsl_matrix *X;
    gsl_matrix *Y;
    double teststat;
} mv_mat;

typedef struct GroupMatrix{
    gsl_matrix *matrix;
} GrpMat;

// ManyGlm related
typedef struct RegressionMethod{
    // regression methods
    uintptr_t model;
    uintptr_t varStab;
    uintptr_t estiMethod;
    double tol;
    uintptr_t n; // used in binomial regression
} reg_Method;

// ManyLM related
// anova.cpp
class AnovaTest {
public: mv_Method *mmRef;
	gsl_matrix *Yref;
	gsl_matrix *Xref;
	gsl_matrix *inRef;
	uintptr_t nSamp;

	double *multstat;
	double *Pmultstat;
	gsl_matrix *statj;
	gsl_matrix *Pstatj;
	uintptr_t *dfDiff;
	gsl_matrix *bootID;	 

       // Methods
       AnovaTest(mv_Method *, gsl_matrix *, gsl_matrix *, gsl_matrix *isXvarIn);
       virtual ~AnovaTest();
       intptr_t resampTest(void); 
       void releaseTest(void);
       void display(void);

private: mv_mat *Hats;
         gsl_permutation **sortid;
	 gsl_vector *bStatj;
	 double bMultStat;
	 uintptr_t nModels, nRows, nVars, nParam;

         // Methods
	 intptr_t anovacase(gsl_matrix *bY, gsl_matrix *bX);
         intptr_t anovaresi(gsl_matrix *bY, const uintptr_t p);
};

// summary.cpp
class Summary
{
public: mv_Method *mmRef;
	gsl_matrix *Yref;
	gsl_matrix *Xref;
	uintptr_t nSamp;

        double R2;	
	double *multstat;
	double *Pmultstat;
	gsl_matrix *unitstat;
	gsl_matrix *Punitstat;
	gsl_matrix *bootID;	 

       // Methods
       Summary(mv_Method *, gsl_matrix *, gsl_matrix *);
       virtual ~Summary();
       intptr_t resampTest(void); 
       void releaseSummary(void);
       void display(void);

private: mv_mat *Hats;
	 gsl_permutation **sortid;
	 uintptr_t nRows, nVars, nParam;
	 double *bMultStat;
	 gsl_matrix *bUnitStat;

	 // Methods
         intptr_t calcR2(void);
	 intptr_t smrycase(gsl_matrix *bY, gsl_matrix *bX);
         intptr_t smryresi(gsl_matrix *bY);

};

// glm base class
class glm
{
   public:
           glm(const reg_Method *mm);
	   virtual ~glm();	   
	   void initialGlm(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B);
	   void releaseGlm(void);
	   virtual intptr_t regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B)=0;
	   virtual intptr_t EstIRLS(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B, double *a)=0;
	   intptr_t copyGlm(glm *src);
           void display(void); 	  

           // input arguments
           const reg_Method *mmRef;
           gsl_matrix *Yref;
	   gsl_matrix *Xref;
	   gsl_matrix *Oref;
	   // return properties
	   gsl_matrix *Beta;
           gsl_matrix *varBeta; // variance of Beta Hat
	   gsl_matrix *Mu;
	   gsl_matrix *Eta;
	   gsl_matrix *Res;
	   gsl_matrix *Var;
	   gsl_matrix *wHalf;
	   gsl_matrix *sqrt1_Hii;
           gsl_matrix *PitRes;

           uintptr_t n; // used in binomial and logistic regression
           uintptr_t rdf;
	   double *phi, *ll, *dev, *aic;
	   uintptr_t *iterconv;  
           uintptr_t maxiter;
	   double eps, maxtol;
           uintptr_t nRows, nVars, nParams;
//   private: 
  	   // abstract 	
	   virtual double link(double) const=0;	   
	   virtual double invLink(double) const=0;	   
	   virtual double rcpLinkDash(double) const=0;	   
	   virtual double weifunc(double, double) const=0;   
	   virtual double varfunc(double, double) const=0;	   
	   virtual double llfunc(double, double, double) const=0;	   
	   virtual double devfunc(double, double, double) const=0;	   
	   virtual double pdf(double, double, double) const=0;	   
	   virtual double cdf(double, double, double) const=0;	   
	   virtual uintptr_t cdfinv(double, double, double) const=0;	   
           virtual uintptr_t genRandist(double, double) const=0;
};

// poisson regression
class PoissonGlm : public glm 
{
    public: // public functions
           PoissonGlm(const reg_Method *);
	   virtual ~PoissonGlm();
	   virtual intptr_t regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B) 
	        { return EstIRLS ( Y, X, O, B, NULL); }
	   intptr_t EstIRLS( gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, double * );
	   intptr_t betaEst( uintptr_t id, uintptr_t iter, double *tol, double a );
	   double getDisper( uintptr_t id ) const;
           intptr_t update(gsl_vector *bj, uintptr_t id);
           intptr_t predict(gsl_vector_view bj, gsl_vector *coef_old, uintptr_t id, double a);


//    private: 
           // Log-link and property functions
           double link(double mui) const 
                { return log(mui); } 
           double invLink(double etai) const 
	        { return exp(etai); }
	   double rcpLinkDash(double mui) const
                { return mui; }
	   double weifunc(double mui, double a) const
	        { return mui; }
	   double varfunc(double mui, double a) const
	        { return mui; }
	   double llfunc(double yi, double mui, double a) const
	        { if (yi==0) return 2*(-mui); else
		  return 2*(yi*log(mui) - mui - gsl_sf_lngamma(yi+1) ); }
	   double devfunc(double yi, double mui, double a) const
	        { return 2*(yi*log(GSL_MAX(yi, eps)/mui)-yi+mui); }
	   double pdf(double yi, double mui, double a) const
                { return Rf_dpois(yi, mui, FALSE); }
	   double cdf(double yi, double mui, double a) const
                { return Rf_ppois(yi, mui, TRUE, FALSE); }
           uintptr_t cdfinv(double u, double mui, double a) const
                { return (uintptr_t)Rf_qpois(u, mui, TRUE, FALSE); }
           uintptr_t genRandist(double mui, double a) const
                { return Rf_rpois(mui); }

};


// Binomial regression yi ~ BIN(n, pi), mui=n*pi
class BinGlm : public PoissonGlm
{
    public: // public functions
           BinGlm(const reg_Method *);
           virtual ~BinGlm();
//    private: // logit link and property functions          
           double link(double mui) const // pi=mui/n
               { return log(mui/(n-mui)); }
           double invLink(double ei) const
                { return n*exp(ei)/(1+exp(ei)); }
           double rcpLinkDash(double mui) const
                { return (mui/n)*(1-mui/n); }
           double weifunc(double mui, double a) const
                { return rcpLinkDash(mui)/n; }
           double varfunc(double mui, double a) const
                { return mui*(1-mui/n); } // n*pi*(1-pi)
           double llfunc(double yi, double mui, double a) const
                { return 2*((yi>0)?yi*log(mui/n):0+(yi<n)?(n-yi)*log(1-mui/n):0+(n>1)?(gsl_sf_lngamma(n+1)-gsl_sf_lngamma(yi+1)-gsl_sf_lngamma(n-yi+1)):0); }
           double devfunc(double yi, double mui, double a) const
                { return 2*((yi>0)?(yi*log(yi/mui)):0+(yi<n)?((n-yi)*log((n-yi)/(n-mui))):0);}
	   double pdf(double yi, double mui, double a) const
//	        { if (n==1) return (yi<1)?(1-mui):mui;
                { return Rf_dbinom(yi, n, mui/n, FALSE); }
           double cdf(double yi, double mui, double a) const
//                { if (n==1) return (yi<1)?(1-mui):1;
                { return Rf_pbinom(yi, n, mui/n, TRUE, FALSE); }
           uintptr_t cdfinv(double ui, double mui, double a) const
//                { if (n==1) return (ui<1)?0:1;
                { return (uintptr_t) Rf_qbinom(ui, n, mui/n, TRUE, FALSE); }
           uintptr_t genRandist(double mui, double a) const
                { return (uintptr_t) Rf_rbinom(n, mui/n); }

};

// negative binomial regression
class NBinGlm : public PoissonGlm // Y~NB(n, p)
{
    public:
           NBinGlm(const reg_Method *mm);
	   virtual ~NBinGlm();
	   virtual intptr_t regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B)
	      { return nbinfit ( Y, X, O, B); }
	   intptr_t nbinfit(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);
//    private:
	   double weifunc(double mui, double a) const
	        { return (mui/(1+a*mui)); }
	   double varfunc(double mui, double a) const
	        { return (mui*(1+a*mui)); }
	   double llfunc(double yi, double mui, double a) const;
	   double devfunc(double yi, double mui, double a) const
                { return 2*((yi==0)?0:yi*log(yi/mui)-(yi+1/GSL_MAX(a,eps))*log((1+yi*a)/(1+mui*a))); }
           // mu=y*p(y), y=k, r=1/phi, p=1-(1+phi*mu)^-1 (the wiki definition)
           // mu=y*p(y), y=k, n=1/phi, p=(1+phi*mu)^-1 (the GSL definition)
	   double pdf(double yi, double mui, double a) const
                { if (a==0) return Rf_dpois(yi, mui, FALSE);
                  else return Rf_dnbinom(yi, 1/a, 1/(1+a*mui), FALSE); }
	   double cdf(double yi, double mui, double a) const
                { if (a==0) return Rf_ppois(yi, mui, TRUE, FALSE);
                  else return Rf_pnbinom(yi, 1/a, 1/(1+a*mui), TRUE, FALSE); }
           uintptr_t cdfinv(double ui, double mui, double a) const
                { if (a==0) return (uintptr_t) Rf_qpois(ui, mui, TRUE, FALSE);
                  else return (uintptr_t)Rf_qnbinom(ui, 1/a, 1/(1+a*mui), TRUE, FALSE); }
           uintptr_t genRandist(double mui, double a) const
                { if (a==0) return (uintptr_t) Rf_rpois(mui);
                  else return (uintptr_t) Rf_rnbinom(1/a, 1/(1+a*mui)); }

	   intptr_t getfAfAdash (double a, uintptr_t id, double *fA, double *fAdash);
};


// base test class
class GlmTest
{
    public: 
            const mv_Method *tm;
            glm *fit;	   
	    gsl_matrix *Xin;
	    gsl_matrix *smryStat, *Psmry;
	    gsl_matrix *anovaStat, *Panova;
	    gsl_matrix *bootID;
            uintptr_t nSamp;	    
            double *aic;
            uintptr_t *dfDiff;

            // methods
            GlmTest(const mv_Method *tm);
            virtual ~GlmTest(void);
            void releaseTest(void);	

	    intptr_t summary(glm *);
	    void displaySmry(glm *);

	    intptr_t anova(glm *, gsl_matrix *);
            void displayAnova(void);
            
    private:
	    intptr_t getBootID(void);

//	    intptr_t geeCalc(glm *PtrAlt, glm *PtrNull, gsl_matrix *);
	    intptr_t GeeWald(glm *, gsl_matrix *, gsl_vector *, double lambda);
	    intptr_t GeeScore(gsl_matrix *, glm *, gsl_vector *, double lambda);
	    intptr_t GeeLR(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat);

//          intptr_t resampSmryCase(glm *, gsl_matrix *, GrpMat *, GrpMat *, uintptr_t i ); // summary
            intptr_t resampSmryCase(glm *, gsl_matrix *, GrpMat *, gsl_matrix *, uintptr_t i ); // summary
	    intptr_t resampAnovaCase(glm *, gsl_matrix *, gsl_matrix *, gsl_matrix *, uintptr_t i);
	    intptr_t resampNonCase(glm *, gsl_matrix *, uintptr_t i);
	    intptr_t setMonteCarlo(glm *model, gsl_matrix *, gsl_matrix *);

	    // the following used in resampling
	    uintptr_t nModels;
            gsl_rng *rnd;
            uintptr_t *permid;   // only useful in permutation test
            double lambda, eps;    // intermediate shrinkage parameter
                
            // the following are used in geeCalc
	    gsl_matrix *L;  // only useful in Wald test
	    gsl_matrix *Rlambda; 
	    gsl_matrix *Wj;
	    gsl_matrix *XBeta, *Sigma; // used in monte carlo simulation
	    //double *mr, *sr; // mean and variance of model residuals
//	    gsl_vector *mr;

	    GrpMat *GrpXs;  // group X0
	    GrpMat *GrpOs;  // group offset

};

// io.cpp - input/output functions
intptr_t vector_filesize(FILE *f);
void matrix_filesize(FILE *f, intptr_t * row, intptr_t * col);
gsl_matrix * load_m(const char * file);
gsl_vector * load_v(const char * file);
void displaymatrix(gsl_matrix * m, const char * name);
void displayvector(gsl_vector * v, const char * name);

// calctest.c - utility functions
double calcDet(gsl_matrix *SS);
intptr_t is_sym_matrix(const gsl_matrix *mat);
intptr_t subX(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
intptr_t subX1(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
intptr_t subX2(gsl_matrix *X, uintptr_t id, gsl_matrix *Xi);
intptr_t subXrow(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
intptr_t subXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
intptr_t addXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
intptr_t subXrow1(gsl_matrix *X, gsl_vector *ref0, gsl_vector *ref1, gsl_matrix *Xi);
intptr_t GetR(gsl_matrix *Res, uintptr_t corr, double lambda, gsl_matrix *R); 
intptr_t subtractMean(gsl_matrix *dat);
// calctest.c - manylm related functions
intptr_t testStatCalc(mv_mat *H0, mv_mat *H1, mv_Method *mmRef, const uintptr_t ifcalcH1det, double *stat, gsl_vector *statj);
intptr_t calcSS(gsl_matrix *Y, mv_mat *Hat, mv_Method *mmRef, const uintptr_t ifcalcHat, const uintptr_t ifcalcCoef, const uintptr_t ifcalcSS);
intptr_t calcAdjustP(const uintptr_t punit, const uintptr_t nVars, double *bj, double *sj, double *pj, gsl_permutation *sortid);
intptr_t reinforceP(double *p, uintptr_t nVars, gsl_permutation *sortid);
intptr_t getHat(gsl_matrix *X, gsl_matrix *W, gsl_matrix *Hat);
intptr_t invLSQ(gsl_matrix *A, gsl_vector *b, gsl_vector *x);
intptr_t rcalc(gsl_matrix *Res, double, uintptr_t, gsl_matrix *SS);

// simutility.cpp - functions used in simulation tests 
intptr_t GetMean(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *Mu);
intptr_t GetPdstbtion(double *p, uintptr_t nVars, uintptr_t *isH0var, uintptr_t *cnt, uintptr_t *cntfwe);
//intptr_t GetCov (gsl_matrix *Mu, gsl_matrix *Y, uintptr_t AR1MAT, gsl_matrix *Sigma);
//intptr_t GetMeanCov(gsl_matrix *X, gsl_matrix *Y, mv_Method *mm, uintptr_t AR1MAT, gsl_matrix *Mu, gsl_matrix *Sigma);

// rnd.c - functions to generate random numbers from multivariate (normal) distributions
// MVN random number generator
intptr_t rmvnorm(const gsl_rng *, const uintptr_t, const gsl_matrix *, gsl_vector *);
// MVN with positive-semi definite covariance matrix
intptr_t semirmvnorm(const gsl_rng *, const uintptr_t, const gsl_matrix *, gsl_vector *);

#endif
