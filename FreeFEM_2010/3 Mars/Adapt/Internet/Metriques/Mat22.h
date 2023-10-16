/*
 *  Mat22.h
 *  
 *All routines for 2x2 symmetric or non symmetric matrices.
 *  Created by Jean-Marie Mirebeau on 05/06/09.
 *  Copyright 2009 UPMC. All rights reserved.
 *
 */
void SquareSym2(const double metric[3], double squareMetric[3]);
int SqrtSym2(const double metric[3], double sqrtMetric[3]);
int AbsSym2(const double metric[3], double absMetric[3]);
int SqrtSym2Safe(const double metric[3], double sqrtMetric[3]);
int InvSym2(const double metric[3], double invMetric[3]);
void BoundExcentricity(const double metric[3], double boundedMetric[3], double maxCond);
void SquareToSym2(const double A[4], double metric[3]);
int InvSquare2(const double A[2][2], double invA[2][2]);
double CondSym2(const double metric[3]);
void ProdSquare2(const double A[2][2], const double v[2], double w[2]);
	
//square of a 2x2 symmetric matrix.
void SquareSym2(const double metric[3], double squareMetric[3]){
	squareMetric[0] = square(metric[0])+square(metric[1]);
	squareMetric[1] = metric[1]*(metric[0]+metric[2]);
	squareMetric[2] = square(metric[2])+square(metric[1]);
};

//square root of a 2x2 symmetric matrix, with error message if non positive.
//in that case, the sqrt of the absolute value is sent.
int SqrtSym2(const double metric[3], double sqrtMetric[3]){
	double Det = metric[0]*metric[2]-square(metric[1]);
	double HTr = (metric[0]+metric[2])/2.;
	
	if(Det<0 or HTr<0){	return 1;}; //On pourrait aussi renvoyer la racine de la vlaeur absolue.
	
	double sqrtDet=sqrt(Det);
	double delta = sqrt(square(HTr)-Det);
	double SumSqrtLambda = sqrt(HTr+delta)+sqrt(HTr-delta);
	
	//la racine carrée de metric est sqrtMetric = (metric+sqrtDet* Id)/SumSqrtLambda.
	if (SumSqrtLambda!=0) {
		sqrtMetric[0]= (metric[0] +sqrtDet)/SumSqrtLambda;
		sqrtMetric[1]=  metric[1]/SumSqrtLambda;
		sqrtMetric[2]= (metric[2]+sqrtDet)/SumSqrtLambda;		
	} else {
		sqrtMetric[0]=0.;
		sqrtMetric[1]=0.;
		sqrtMetric[2]=0.;
	};
	
	return 0;
};

//absolute value of a metric.
//returns 0 if the metric was non-negative, and -1 otherwise
int AbsSym2(const double metric[3], double absMetric[3]){
	if(metric[0]*metric[2]-square(metric[1]) >= 0 and metric[0]+metric[2]>=0){
		absMetric[0] = metric[0];
		absMetric[1] = metric[1];
		absMetric[2] = metric[2];
		return 0;
	};
	
	double squareMetric[3];
	SquareSym2(metric, squareMetric);
	SqrtSym2(squareMetric, absMetric);
	return -1;
};

//square root of the absolute value of a 2x2 matrix.
int SqrtSym2Safe(const double metric[3], double sqrtMetric[3]){
	double absMetric[3];
	int sign = AbsSym2(metric, absMetric);
	SqrtSym2(absMetric, sqrtMetric);
	return sign;
};


int InvSym2(const double metric[3], double invMetric[3]){
	double det = metric[0]*metric[2]-square(metric[1]);
	if(det==0.){invMetric[0]=1.; invMetric[1]=0.; invMetric[2]=1.; return 1;};
	invMetric[0] = metric[2]/det;
	invMetric[1] = -metric[1]/det;
	invMetric[2] = metric[0]/det;
	return 0;
};

int InvSquare2(const double A[2][2], double invA[2][2]){
	double det = A[0][0]*A[1][1]-A[1][0]*A[0][1];
	if(det==0.){invA[0][0]==1; invA[1][0]=0.; invA[0][1]=0.; invA[1][1]=1.; return 1;};
	invA[0][0]=  A[1][1]/det;
	invA[1][0]= -A[1][0]/det;
	invA[0][1]= -A[0][1]/det;
	invA[1][1]=  A[0][0]/det;
	return 0;
};

//int BoundExcentricity...
/*La routine suivante sert à borner l'excentricité d'une matrice.
 plus précisément, étant donnée une matrice symétrique définie positive, on remplace la plus petite vap
 par un truc plus grand pour contenir le condtionnement.
 */
void BoundExcentricity(const double metric[3], double boundedMetric[3], double maxCond){
	//boundedMetric[0]=metric[0]; boundedMetric[1]=metric[1]; boundedMetric[1]=metric[1];//pour pouvoir quitter en cas de pb.
	//if(maxCond<=0){return 1;};
	
	double det = metric[0]*metric[2]-square(metric[1]);
	double hTr = (metric[0]+metric[2])/2;
	double hDiff = sqrt(square(hTr)-det); //en théorie c'est positif, si ma matrice est positive.
	double lambdam = hTr-hDiff;
	double lambdaM = hTr+hDiff;
	double lambdaN = max(lambdam,lambdaM/max(1.,maxCond)); //on borne le rapport des vap.
	
	if(lambdam == lambdaM){boundedMetric[0]=lambdaM; boundedMetric[1]=0; boundedMetric[2]=lambdaM; return;};
	
	double a=(lambdaM-lambdaN)/(lambdaM-lambdam);
	double b=lambdaM*(lambdaN-lambdam)/(lambdaM-lambdam);
	
	boundedMetric[0] = a*metric[0]+b;
	boundedMetric[1] = a*metric[1];
	boundedMetric[2] = a*metric[2]+b;
	return;
};

void SquareToSym2(const double A[2][2], double metric[3]){//calcule transpose(A).A
	metric[0] = square(A[0][0]) + square(A[1][0]);
	metric[1] = A[0][0]* A[0][1] + A[1][0]*A[1][1];
	metric[2] = square(A[0][1]) + square(A[1][1]);
};


void SquareToSym2_AAT(const double A[2][2], double metric[3]){//calcule A.transpose(A)
	metric[0] = square(A[0][0]) + square(A[0][1]);
	metric[1] = A[0][0]* A[1][0] + A[0][1]*A[1][1];
	metric[2] = square(A[1][0]) + square(A[1][1]);
};


double CondSym2(const double metric[3]){
	double hSum = (metric[0]+metric[2])/2.;
	double hDiff = sqrt(square(hSum) - (metric[0]*metric[2]-square(metric[1])) );
	if(hSum<=hDiff){return TGV;};
	return (hSum+hDiff)/(hSum-hDiff);
};

void ProdSquare2(const double A[2][2], const double v[2], double w[2]){
	w[0] = A[0][0]*v[0] +A[0][1]*v[1];
	w[1] = A[1][0]*v[0] +A[1][1]*v[1];
};