/*
 *  BasicMath.h
 *  AllMetricsToCL
 *
 *  Created by Jean-Marie Mirebeau on 10/06/09.
 *  Copyright 2009 UPMC. All rights reserved.
 *
 */

#define TGV 1000000.
#define TPV 1/1000000.

inline double square(double a){return a*a;};
inline double cube(double a){return a*a*a;};
inline double sign(double a){return a>0? 1: (a<0? -1 :0);};

double argmin(double table[], int len){
	double min=table[0]; int pos=0;
	for(int i=1; i<len; i++){if(table[i]<min){min = table[i]; pos=i;};};
	return pos;
};

double argmax(double table[], int len){
	double max=table[0]; int pos=0;
	for(int i=1; i<len; i++){if(table[i]>max){max = table[i]; pos=i;};};
	return pos;
};
