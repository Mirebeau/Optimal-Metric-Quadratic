/*
 *  
 *  New metrics for FreeFem++
 *
 *  Created by Jean-Marie Mirebeau on 05/06/09.
 *  Copyright 2009 UPMC. All rights reserved.
 *
 */

#include "config.h"
#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "AFunction_ext.hpp"
#include <cstdlib>

/*******************Les fonctions contenues dans ce document************/
inline double square(double u){return u*u;};

double DetSym(const double S[3]);
double NormSym(const double S[3]);
void MultSym(double S[3], double lambda);
void MultSymDet(double S[3], double alpha);
void AffSym(double S[3], double a, double b);
void SquareSym(double S[3]);
void EigenSym(const double S[3], double E[2]);
void SqrtSym(double S[3]);
void AbsSym(double S[3]);
void MaxSym(double S[3], double lambda);
void BoundExcentricity(double S[3], double maxCond);
void CopySym(const double SSource[3], double SDest[3]);
double DiscCube(const double D[4]);
double NormCube(const double D[4]);
double RatioCube(const double D[4]);

/*******************les métriques***************/

//métrique pour la norme H^1, les éléments P_1, robuste aux angles obtus.
void MetricObP1(const double H[3], double M[3], double maxCond=-1.){
	CopySym(H,M);
	AbsSym(M);
	MultSym(M, NormSym(M));
	BoundExcentricity(M,  maxCond);	
	MultSymDet(M,-1/4.);
};

void MetricObP2(const double D[4], double M[3], double maxCond=-1.){
	M[0] = square(D[0])+2*square(D[1])+square(D[2]);
	M[1] = D[0]*D[1]+2*D[1]*D[2]+D[2]*D[3];
	M[2] = square(D[1])+2*square(D[2])+square(D[3]);
	SqrtSym(M);
	MaxSym(M,RatioCube(D));
	BoundExcentricity(M,  maxCond);	
	MultSymDet(M,-1/6.);
};

/*****************La manipulation des matrices 2x2*****************/
double DetSym(const double S[3]) {return S[0]*S[2]-square(S[1]);};

double NormSym(const double S[3]) {return fabs(S[0]+S[2])/2. + sqrt( square(S[0]-S[2])/4 + square(S[1]));};
/*{
	double det = metric[0]*metric[2]-square(metric[1]);
	double hTr = (metric[0]+metric[2])/2;
	double hDiff = sqrt(square(hTr)-det); 
	return fabs(hTr)+hDiff;
};*/

void MultSym(double S[3], double lambda){
	S[0] = lambda*S[0];
	S[1] = lambda*S[1];
	S[2] = lambda*S[2];
};

void MultSymDet(double S[3], double alpha){
	double det = DetSym(S);
	if(det>0 || alpha>=0) MultSym(S, pow(det, alpha));
};


//image d'une matrice symétrique par une fonction affine.
void AffSym(double S[3], double a, double b){
	S[0] = a*S[0]+b;
	S[1] = a*S[1];
	S[2] = a*S[2]+b;
};

void SquareSym(double S[3]){
	double SS[3];
	SS[0] = square(S[0])+square(S[1]);
	SS[1] = S[1]*(S[0]+S[2]);
	SS[2] = square(S[1])+square(S[2]);
	CopySym(SS,S);
};

//E[0] : petite vap, E[1] : grande vap de S
void EigenSym(const double S[3], double E[2]){
	double hDiff = sqrt( square(S[0]-S[2])/4 + square(S[1]));
	double hSum = (S[0]+S[2])/2;
	E[0] = hSum-hDiff;
	E[1] = hSum+hDiff;
	return;
};

void SqrtSym(double S[3]){
	double E[2];
	EigenSym(S, E);
	if(E[0]<0 || E[1]<=0) return;
	double a = 1/(sqrt(E[0])+sqrt(E[1]));
	double b = a*sqrt(E[0]*E[1]);
	AffSym(S,a,b);
};

void AbsSym(double S[3]){
	SquareSym(S);
	SqrtSym(S);
};

void MaxSym(double S[3], double lambda){
	double E[2];
	EigenSym(S,E);
	if(lambda <= E[0]) return;
	if(E[1]	  <= lambda) {S[0] = lambda; S[1]=0; S[2] = lambda; return;};
	AffSym(S, (E[1]-lambda)/(E[1]-E[0]), E[1] * (lambda - E[0])/(E[1] - E[0]));
};
	
void BoundExcentricity(double S[3], double maxCond){
	if(maxCond>=1)	MaxSym(S, NormSym(S)/square(maxCond)); 	
};


void CopySym(const double SSource[3], double SDest[3]){
	SDest[0] = SSource[0]; SDest[1] = SSource[1]; SDest[2] = SSource[2];
};
		   
double DiscCube(const double D[4]) {return 27*(4*(D[0]*D[2]-square(D[1]))*(D[1]*D[3]-square(D[2])) - square(D[0]*D[3]-D[1]*D[2]));};
double NormCube(const double D[4]) {return sqrt(square(D[0])+3*square(D[1])+3*square(D[2])+square(D[3]));};
double RatioCube(const double D[4]){
	double n=NormCube(D);
	if(n>0) return pow(max(0.,-DiscCube(D))/n, 1/3.);
	return 0.;
};
		   
		   
/*****************L'interface avec FreeFem*************/

/*********************** P1H1 **********************/
double MetricP1xx(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricObP1(deriv,metric);
	return metric[0];
};
double MetricP1xy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricObP1(deriv,metric);
	return metric[1];
};
double MetricP1yy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricObP1(deriv,metric);
	return metric[2];
};

//P1H1, avec conditionnement
double MetricP1Condxx(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricObP1(deriv,metric, cond);
	return metric[0];
};
double MetricP1Condxy(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricObP1(deriv,metric, cond);
	return metric[1];
};
double MetricP1Condyy(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricObP1(deriv,metric, cond);
	return metric[2];
};

/*************** P2H1 ***************/
double MetricP2xx(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricObP2(deriv,metric);
	return metric[0];
};
double MetricP2xy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricObP2(deriv,metric);
	return metric[1];
};
double MetricP2yy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricObP2(deriv,metric);
	return metric[2];
};

//P2H1 avec conditionnement
double MetricP2Condxx(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricObP2(deriv,metric,cond);
	return metric[0];
};
double MetricP2Condxy(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricObP2(deriv,metric,cond);
	return metric[1];
};
double MetricP2Condyy(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricObP2(deriv,metric,cond);
	return metric[2];
};

/************************ Classe nécessaire pour FreeFem++ *********************/

class Init { public:
	Init();
};
Init init;

Init::Init(){
	
	/****P1H1****/
	Global.Add("MetricP1xx","(",new OneOperator3_<double>(MetricP1xx));
	Global.Add("MetricP1xy","(",new OneOperator3_<double>(MetricP1xy));
	Global.Add("MetricP1yy","(",new OneOperator3_<double>(MetricP1yy));
	
	Global.Add("MetricP1Condxx","(",new OneOperator4_<double>(MetricP1Condxx));
	Global.Add("MetricP1Condxy","(",new OneOperator4_<double>(MetricP1Condxy));
	Global.Add("MetricP1Condyy","(",new OneOperator4_<double>(MetricP1Condyy));
		
	/****P2H1****/
	Global.Add("MetricP2xx","(",new OneOperator4_<double>(MetricP2xx));
	Global.Add("MetricP2xy","(",new OneOperator4_<double>(MetricP2xy));
	Global.Add("MetricP2yy","(",new OneOperator4_<double>(MetricP2yy));
	
	Global.Add("MetricP2Condxx","(",new OneOperator5_<double>(MetricP2Condxx));
	Global.Add("MetricP2Condxy","(",new OneOperator5_<double>(MetricP2Condxy));
	Global.Add("MetricP2Condyy","(",new OneOperator5_<double>(MetricP2Condyy));	
};
