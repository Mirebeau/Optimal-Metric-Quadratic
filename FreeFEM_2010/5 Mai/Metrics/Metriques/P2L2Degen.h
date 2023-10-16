/*
 *  P2L2Degen.h
 *  
 *
 *  Created by Jean-Marie Mirebeau on 03/09/09.
 *  Copyright 2009 UPMC. All rights reserved.
 *
 */
const double CTripleRoot = 1e-2;
double distP1(double u, double v);
//input : a (close to) degenerate polynomial deriv, (close to) conjugated to x^2 y or x^3
//output : the optimal metric with prescribed conditionning 
void MetricP2L2Degen(const double deriv[4], double metric[3], double maxCond=-1.){
	double roots[3];
	double phi[2][2]; double invphi[2][2];
	double normPi;
	double D[2];
	double sign = ThirdDegreeRealDegen(deriv[0], 3*deriv[1], 3*deriv[2], deriv[3], roots);
	
	//metric[0]=0; metric[1] = 1; metric[2] = 2; return;
	//metric[0] = deriv[0]; metric[1] = deriv[1]; metric[2] = deriv[3]; return;
	//metric[0] = roots[0]; metric[1] = roots[1]; metric[2] = roots[2]; return;
	
	if(sign<0){roots[2]=roots[1];}; //on oublie les parties complexes
	
	double dist[3]; 
	dist[0]=distP1(roots[0],roots[1]);  
	dist[1]=distP1(roots[1],roots[2]); 
	dist[2]=distP1(roots[2],roots[0]);
	
	int amin = argmin(dist,3);
	double dmax = max(dist[0],max(dist[1],dist[2]));
	
	//metric[0] = dmax;
	//return;
	
	//metric[0] = roots[0]; metric[1] = roots[1]; metric[2] = roots[2]; return;
	
	//metric[0] = roots[amin]; metric[1] = roots[(amin+2)%3]; return;
	//metric[0] = dmax; metric[1] = distP1(2,2); metric[2] = dist[2]; return;
	
	if(dmax < CTripleRoot){ //trois racines identiques. Correct ?
		Rotation(atan(roots[0]), phi);
		//gérer maxCond<1???
		D[0]=1; D[1] = 1/maxCond;
		SquareDiagToSym2(phi, D, metric);
		//la norme invariante par rotation
		normPi = square(deriv[0])+3*square(deriv[1])+3*square(deriv[2])+square(deriv[3]);
		normPi = sqrt3(normPi);
		metric[0]*= normPi; 
		metric[1]*= normPi; 
		metric[2]*= normPi; 
		return;
	};
	
	double DRoot = roots[amin]; //double and simple root
	double SRoot = roots[(amin+2)%3];
	
	//metric[0] = DRoot; metric[1] = SRoot; return;
	
	BiRotation(M_PI/2.-atan(SRoot),M_PI/2.-atan(DRoot),invphi); //pour la métrique
	InvSquare2(invphi, phi); //pi o invphi = x^2 y
	
	//metric[0] = phi[0][0]; metric[1]=phi[0][1]; metric[2] = phi[1][0]; return;
	//coef dominant de pi o phi= lambda x^2 y à récupérer.
	//pour cela on teste sur le vecteur (1,1)
	double v[2] ={phi[0][0]+phi[0][1], phi[1][0]+phi[1][1]};
	normPi = deriv[0]*cube(v[0]) +3*deriv[1]*square(v[0])*v[1] +3*deriv[2]*v[0]*square(v[1]) +deriv[3]*cube(v[1]);
	
	//metric[0] = v[0]; metric[1] = v[1]; metric[2] = normPi; return;
	
	//Solve for conditionning
	double trM1 = square(invphi[0][0])+square(invphi[0][1]);
	double trM2 = square(invphi[1][1])+square(invphi[1][0]); trM2*=4./27.;
	double detM = square(invphi[0][0]*invphi[1][1]-invphi[1][0]*invphi[0][1])*4./27.;
	double ratioCond = square(sqrt(maxCond)+1/sqrt(maxCond));
	
	//metric[0]=trM1; metric[1]=trM2; metric[2]=detM; return; 
	
	//eqn (trM1*lambda+trM2)^2 = ratioCond*lambda*detM
	//i.e. lambda^2*trM1^2 + lambda*(2*trM1*trM2 - ratioCond*detM)+trM2^2 == 0;
	//normalement deux racines positives, et on prend la plus petite.
	//sinon deux racines complexes conjuguées. On prend l'ellipsoide critique 
	//puis on applique la routine pour borner l'excentricité.
	//eqn (trM1*trM2-ratioCond*detM/2)^2 - square(trM1)square(trM2)
	// se factorise. racine double assoc : lambdaMin
	double ratioCondMin = 4*trM1*trM2/detM;
	double lambdaMin = trM2/trM1; double lambda;
	
	if(ratioCond>ratioCondMin){
		double lambdaList[2];
		SecondDegreeReal(square(trM1), 2*trM1*trM2-ratioCond*detM, square(trM2),lambdaList);
		//metric[0] = lambdaList[0]; metric[1]=lambdaList[1];return;
		lambda = lambdaList[0]>lambdaList[1]? sqrt3(lambdaList[0]) : sqrt3(lambdaList[1]);
		D[0] = lambda; D[1] = 4/(27*square(lambda));
		SquareDiagToSym2(invphi, D, metric);
	} else { 
		lambda = sqrt3(lambdaMin);
		D[0] = lambda; D[1] = 4/(27*square(lambda));
		double metricLong[3];
		SquareDiagToSym2(invphi, D, metricLong);
		BoundExcentricity(metricLong, metric, maxCond);
	};
	
	//metric[0] = lambda; return;
	
	normPi = pow(fabs(normPi),2/3.);
	metric[0]*=normPi;
	metric[1]*=normPi;
	metric[2]*=normPi;
	return;
	
};


//on est obligé de faire une nouvelle routine po
