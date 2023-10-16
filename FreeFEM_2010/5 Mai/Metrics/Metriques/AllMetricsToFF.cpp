/*
 *  AllMetrics.h
 *  
 *  This file adds all developped metrics to FreeFem++
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
#include <complex>

/**** Des routines nécessaires ****/
#include "BasicMath.h"
#include "Algebra.h" //utile seulement pour P2L2Cond
#include "Mat22.h"

/**** Les métriques ****/
#include "P2L2Cond.h"
#include "P1L2H1.h"

/*********************** Interfaçage avec FreeFem++ **********************/

/****** On prend les coordonnées des métriques et on les renvoie à FF une par une *******/
// cela fait calculer les métriques trois fois mais tant pis.

/*********************** P1L2 **********************/
double MetricP1L2xx(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1L2(deriv,metric);
	return metric[0];
};

double MetricP1L2xy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1L2(deriv,metric);
	return metric[1];
};
double MetricP1L2yy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1L2(deriv,metric);
	return metric[2];
};

//P1L2 avec Conditionnement
double MetricP1L2Condxx(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1L2(deriv,metric,cond);
	return metric[0];
};

double MetricP1L2Condxy(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1L2(deriv,metric,cond);
	return metric[1];
};
double MetricP1L2Condyy(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1L2(deriv,metric,cond);
	return metric[2];
};

/*********************** P1H1 **********************/
double MetricP1H1xx(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric);
	return metric[0];
};
double MetricP1H1xy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric);
	return metric[1];
};
double MetricP1H1yy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric);
	return metric[2];
};

//P1H1, avec conditionnement
double MetricP1H1Condxx(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric, cond);
	return metric[0];
};
double MetricP1H1Condxy(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric, cond);
	return metric[1];
};
double MetricP1H1Condyy(const double& a, const double& b, const double& c, const double& cond){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric, cond);
	return metric[2];
};

//P1H1, avec erreur due aux grands angles contrôlée.
double MetricP1H1Sxx(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric, -1., 1.);
	return metric[0];
};
double MetricP1H1Sxy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric, -1., 1.);
	return metric[1];
};
double MetricP1H1Syy(const double& a, const double& b, const double& c){
	double deriv[3]; deriv[0]=a; deriv[1]=b; deriv[2]=c;
	double metric[3];
	MetricP1H1(deriv,metric, -1., 1.);
	return metric[2];
};


/**************** P2L2 *****************/
double MetricP2L2xx(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2L2(deriv,metric);
	return metric[0];
};

double MetricP2L2xy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2L2(deriv,metric);
	return metric[1];
};

double MetricP2L2yy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2L2(deriv,metric);
	return metric[2];
};


//P2L2 avec conditionnement
double MetricP2L2Condxx(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2L2(deriv,metric, cond);
	return metric[0];
};

double MetricP2L2Condxy(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2L2(deriv,metric, cond);
	return metric[1];
};

double MetricP2L2Condyy(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2L2(deriv,metric, cond);
	return metric[2];
};

/*************** P2H1 ***************/
double MetricP2H1xx(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric);
	return metric[0];
};
double MetricP2H1xy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric);
	return metric[1];
};
double MetricP2H1yy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric);
	return metric[2];
};

//P2H1 avec conditionnement
double MetricP2H1Condxx(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric,cond);
	return metric[0];
};
double MetricP2H1Condxy(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric,cond);
	return metric[1];
};
double MetricP2H1Condyy(const double& a, const double& b, const double& c, const double& d, const double& cond){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric,cond);
	return metric[2];
};

//P2H1, avec erreur due aux grands angles contrôlée.
double MetricP2H1Sxx(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric, -1., 1.);
	return metric[0];
};
double MetricP2H1Sxy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric, -1., 1.);
	return metric[1];
};
double MetricP2H1Syy(const double& a, const double& b, const double& c, const double& d){
	double deriv[4]; deriv[0]=a; deriv[1]=b; deriv[2]=c; deriv[3]=d;
	double metric[3];
	MetricP2H1(deriv,metric, -1., 1.);
	return metric[2];
};
/************************ Classes nécessaires pour FreeFem++ *********************/

class Init { public:
	Init();
};
Init init;

//le discriminant pour FreeFem++
double Discriminant(const double& a, const double& b, const double& c, const double& d){
	return square(b)*square(c) - 4*a*cube(c) - 4*cube(b)*d + 18*a*b*c*d - 27*square(a)*square(d);
};

Init::Init(){
	/****P1L2****/
	Global.Add("MetricP1L2xx","(",new OneOperator3_<double>(MetricP1L2xx));
	Global.Add("MetricP1L2xy","(",new OneOperator3_<double>(MetricP1L2xy));
	Global.Add("MetricP1L2yy","(",new OneOperator3_<double>(MetricP1L2yy));
	
	Global.Add("MetricP1L2Condxx","(",new OneOperator4_<double>(MetricP1L2Condxx));
	Global.Add("MetricP1L2Condxy","(",new OneOperator4_<double>(MetricP1L2Condxy));
	Global.Add("MetricP1L2Condyy","(",new OneOperator4_<double>(MetricP1L2Condyy));
	
	/****P1H1****/
	Global.Add("MetricP1H1xx","(",new OneOperator3_<double>(MetricP1H1xx));
	Global.Add("MetricP1H1xy","(",new OneOperator3_<double>(MetricP1H1xy));
	Global.Add("MetricP1H1yy","(",new OneOperator3_<double>(MetricP1H1yy));
	
	Global.Add("MetricP1H1Condxx","(",new OneOperator4_<double>(MetricP1H1Condxx));
	Global.Add("MetricP1H1Condxy","(",new OneOperator4_<double>(MetricP1H1Condxy));
	Global.Add("MetricP1H1Condyy","(",new OneOperator4_<double>(MetricP1H1Condyy));

	Global.Add("MetricP1H1Sxx","(",new OneOperator3_<double>(MetricP1H1Sxx));
	Global.Add("MetricP1H1Sxy","(",new OneOperator3_<double>(MetricP1H1Sxy));
	Global.Add("MetricP1H1Syy","(",new OneOperator3_<double>(MetricP1H1Syy));
	
	/****P2L2****/
	Global.Add("Discriminant","(",new OneOperator4_<double>(Discriminant));
	
	Global.Add("MetricP2L2xx","(",new OneOperator4_<double>(MetricP2L2xx));
	Global.Add("MetricP2L2xy","(",new OneOperator4_<double>(MetricP2L2xy));
	Global.Add("MetricP2L2yy","(",new OneOperator4_<double>(MetricP2L2yy));
	
	Global.Add("MetricP2L2Condxx","(",new OneOperator5_<double>(MetricP2L2Condxx));
	Global.Add("MetricP2L2Condxy","(",new OneOperator5_<double>(MetricP2L2Condxy));
	Global.Add("MetricP2L2Condyy","(",new OneOperator5_<double>(MetricP2L2Condyy));
	
	/****P2H1****/
	Global.Add("MetricP2H1xx","(",new OneOperator4_<double>(MetricP2H1xx));
	Global.Add("MetricP2H1xy","(",new OneOperator4_<double>(MetricP2H1xy));
	Global.Add("MetricP2H1yy","(",new OneOperator4_<double>(MetricP2H1yy));
	
	Global.Add("MetricP2H1Condxx","(",new OneOperator5_<double>(MetricP2H1Condxx));
	Global.Add("MetricP2H1Condxy","(",new OneOperator5_<double>(MetricP2H1Condxy));
	Global.Add("MetricP2H1Condyy","(",new OneOperator5_<double>(MetricP2H1Condyy));	

	Global.Add("MetricP2H1Sxx","(",new OneOperator4_<double>(MetricP2H1Sxx));
	Global.Add("MetricP2H1Sxy","(",new OneOperator4_<double>(MetricP2H1Sxy));
	Global.Add("MetricP2H1Syy","(",new OneOperator4_<double>(MetricP2H1Syy));		
}