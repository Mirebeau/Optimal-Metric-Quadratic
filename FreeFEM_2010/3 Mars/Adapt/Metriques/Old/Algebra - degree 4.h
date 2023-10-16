/*
 *  Algebra.h
 *  
 *
 *  Created by Jean-Marie Mirebeau on 04/06/09.
 *  Copyright 2009 UPMC. All rights reserved.
 *
 */


//la partie équations algébriques

const double TGV=1e20;  //très grande valeur
const double TPV=1e-20; //très petite valeur


//résolution de a x +b =0.
complex<double> FirstDegree(complex<double> a, complex<double> b){
	if(abs(a)==0.){ return complex<double>(TGV,0.);};
	return -b/a;
};

//sur les réels
double FirstDegreeReal(double a, double b){
	if(a==0.){ return TGV;};
	return -b/a;
};

//résolution de l'équation de seond degré, le préliminaire indispensable à l'équation de degré 3.
int SecondDegree(complex<double> a, complex<double> b, complex<double> c, complex<double> roots[]){
	
	roots[0] = complex<double>(1.,2.); roots[1] = complex<double>(3.,4.);
	
	if(abs(a)==0.){roots[0] = complex<double>(TGV,0.); roots[1] = FirstDegree(b,c); return 1;};
	complex<double> delta = b*b - (complex<double> (4.,0.))*a*c;
	roots[0] = (-b-sqrt(delta))/((complex<double> (2.,0.))*a);
	roots[1] = (-b+sqrt(delta))/((complex<double> (2.,0.))*a);
	return 0;
};

//faut-il faire un cas particulier pour l'équation réelle? 
//genre renvoyer les racines si elles sont réelles,
//et 0 sinon avec message d'erreur?

int SecondDegreeReal(double a, double b, double c, double roots[]){
	if(a==0.){ roots[0] = TGV; roots[1] = FirstDegreeReal(b, c); return 1;}; //polynôme dégénéré-> code erreur 1.
	double delta = b*b-4*a*c;
	//discriminant négatif-> code erreur -1.
	//On renvoie la partie réelle et imaginaire des deux racines.
	if(delta<0){roots[0] = -b/(2*a); roots[1]== sqrt(-delta); return-1;};  
	
	roots[0] = (-b-sqrt(delta))/(2*a);
	roots[1] = (-b+sqrt(delta))/(2*a);
	return 0;
};

//solutions de l'équation du troisième degré. Les coefficients sont a, b, c, d, la solution doit être renvoyée dans roots. 
//équation a x^3+ b x^2 y+ c x y^2+ d y^3=0
//en fait, a,b,c,d pourraient être simplement réels (seulement ce cas est utile). Mais bon.


int ThirdDegree(complex<double> a, complex<double> b, complex<double> c, complex <double> d, complex<double> roots[]){
	complex<double> rootsSecond[2];
	
	//cas à faire : si a==0. retourner TGV pour l'un, et les solutions de l'équation du second degré pour les autres.
	if(abs(a)==0.){
		roots[2] = (complex<double> (10.,0.));
		SecondDegree(b,c,d,roots);
		return 1;
	};
	
	if(abs(d) == 0.){
		roots[2]==complex<double>(0.,0.);
		SecondDegree(a,b,c,roots);
		return 1;
	};
	
	
	//sinon on normalise, on veut remplacer par l'équation X^3+ p X+q==0, où X=alpha*x+beta.
	complex<double> alpha = pow(a,1/3.); 
	complex<double> beta  = b/(3.*alpha*alpha); //espérons que si a!=0 (numériquement), alors a^2 aussi.
	
	complex<double> p=(c-3.*alpha*beta*beta)/alpha;
	complex<double> q= d-beta*beta*beta-beta*p;
	
	//résolution de l'équation normalisée.
	
	//on calcule les racines de l'équation Y^2+q Y-p^3/27==0, on note r,s leurs racines cubiques.
	//une solution de l'eqn est donnée par r+s, donc une solution de l'eqn originale par (r+s-beta)/alpha.
	
	SecondDegree(complex<double>(1.,0.), q, -p*p*p/complex<double>(27.,0.), rootsSecond);
	
	//à partir d'ici, on a vraiment besoin de nombres complexes.
	
	rootsSecond[0] = pow(rootsSecond[0],1/3.);
	if(abs(rootsSecond[0])!=0){
		rootsSecond[1] = -p/(complex<double>(3.,0.)*rootsSecond[0]);
	} else {
		rootsSecond[1] = pow(rootsSecond[1],1/3.);
	};
	
	complex<double> root1 = rootsSecond[0]+rootsSecond[1];
	root1 = (root1-beta)/alpha;
	
	
	//les autres racines sont obtenues car on connait leur produit  root2 root3 = -d/root1, et leur somme root2+root3 = -b-root1.
	
	complex<double> HSum= (-b/a-root1)/2.;				//demi somme des deux racines restantes.
	complex<double> Hdiff = sqrt(HSum*HSum + d/(a*root1));	//demi différence des deux racines restantes.
	//remarque : on a une division par root1. mais root1 peut être nul, par exemple pour le pol x^3. Dans ce cas, d est nul aussi...
	//le plus simple est d'intercepter ce comportement au début en testant d.
	
	roots[0]=root1;
	roots[1]=HSum+Hdiff;
	roots[2]=HSum-Hdiff;
	
	return 0;
};


/*fonctionnement différent : ne s'applique qu'à des polynômes réels.
Si trois racines réelles, on les renvoie, et on revoie dans roots, et on renvoie 0.
Si une racine réelle et deux complexes, on renvoie la racine réelle, la partie réelle
puis la partie imaginaire des racines complexes, qui sont conjugées.

Notons : faut-il vraiment prendre en compte le coefficient a? (qui peut être normalisé)
 Disons que oui.
*/ 

int ThirdDegreeReal(double a, double b, double c, double d, double roots[]){
	double alpha = -b*b+3*a*c;
	double beta = -3 b*b*b +9*a*b*c - 27*a*a*d;
	double delta = 4*alpha*alpha*alpha +beta*beta;
	
	if(delta>0){ //rappel : delta * disc<0
		complex<double> gamma = pow... //rappel : éviter gamma = 0
		complex<double> gamma2 = alpha/gamma;
		const complex<double> j = complex<double>(-0.5, sqrt(3.)/2);
		const complex<double> j2= complex<double>(-0.5, -sqrt(3.)/2);
		roots[0] = real(-b-gamma2+gamma)/(3*a);
		roots[1] = real(-b  -j*gamma2 +j2*gamma)/(3*a);
		roots[2] = real(-b -j2*gamma2 +j*gamma)/(3*a);
		return 1; //disc positif, trois racines réelles.
	} else {
		double gamma = pow...
		double gamma2 = alpha/gamma;
		roots[0] = (-b-gamma2+gamma)/(3*a);
		roots[1] = (-b +0.5*gamma2-0.5*gamma)/(3*a);
		roots[2] = (gamma2-gamma)/(2*sqrt(3)*a);
		return -1; //disc négatif, 1 réelle et deux complexes.
	};
};

//x^4 +a x^2+b x +c
void FourthDegreeReal3Coef(double a, double b, double c, double roots[]){
	double roots3[3];
	ThirdDegreeReal(-8,4*a,8*c,b^2-4*a*c, roots);
	lambda =roots[0]; //on veut une racine qcq de cette équation.
	
	complex<double> p = sqrt(2*lambda - a);
	complex<double> q=abs(p)>0?-b/(2*p):sqrt(lambda^2-c); //p =0?
	
	//quel est le ormat de retour des racines?
	//ci dessous, ça ne marche certainement pas, car Second Degree renvoie des complexes.
	
	SecondDegree(1, p, q+lambda, roots);
	SecondDegree(1, -p, -q+lambda, roots+2); //on décale de deux pour mettre les racines
};

void FourthDegreeReal4Coef(double b, double c, double d, double e, double roots[]){
	
};

void FourthDegreeReal(double a, double b, double c, double d, double e, double roots[]){
	
	if(a==0.){roots[3]=
	
	double sa = sqrt(fabs(a));
	double ssa = sqrt(sa);
	
	FourthDegreeRealSimple(-3b*b/(8*a*sa)...
	double roots3[3];
	ThirdDegreeReal(-8,4*a,8*c,b^2-4*a*c, roots);
	lambda =roots[0]; //on veut une racine qcq de cette équation.
	
	complex<double> p = sqrt(2*lambda - a);
	complex<double> q=-b/(2*p); //p =0?
	
	SecondDegree(1, p, q+lambda, roots);
	SecondDegree(1, -p, -q+lambda, roots+2);
};



/*	
int ThirdDegreeReal(double a, double b, double c, double d, double roots[]){
	complex<double> rootsSecond[2];
	
	//cas à faire : si a==0. retourner TGV pour l'un, et les solutions de l'équation du second degré pour les autres.
	if(abs(a)==0.){
		roots[2] = (complex<double> (TGV,0.));
		SecondDegree(b,c,d,roots);
		return 1;
	};
	
	if(abs(d) == 0.){
		roots[2]==complex<double>(0.,0.);
		SecondDegree(a,b,c,roots);
		return 1;
	};
	
	
	//sinon on normalise, on veut remplacer par l'équation X^3+ p X+q==0, où X=alpha*x+beta.
	complex<double> alpha = pow(a,1/3.); 
	complex<double> beta  = b/(3.*alpha*alpha); //espérons que si a!=0 (numériquement), alors a^2 aussi.
	
	complex<double> p=(c-3.*alpha*beta*beta)/alpha;
	complex<double> q= d-beta*beta*beta-beta*p;
	
	//résolution de l'équation normalisée.
	
	//on calcule les racines de l'équation Y^2+q Y-p^3/27==0, on note r,s leurs racines cubiques.
	//une solution de l'eqn est donnée par r+s, donc une solution de l'eqn originale par (r+s-beta)/alpha.
	
	SecondDegree(complex<double>(1.,0.), q, -p*p*p/complex<double>(27.,0.), rootsSecond);
	
	//à partir d'ici, on a vraiment besoin de nombres complexes.
	
	rootsSecond[0] = pow(rootsSecond[0],1/3.);
	if(abs(rootsSecond[0])!=0){
		rootsSecond[1] = -p/(complex<double>(3.,0.)*rootsSecond[0]);
	} else {
		rootsSecond[1] = pow(rootsSecond[1],1/3.);
	};
	
	complex<double> root1 = rootsSecond[0]+rootsSecond[1];
	root1 = (root1-beta)/alpha;
	
	
	//les autres racines sont obtenues car on connait leur produit  root2 root3 = -d/root1, et leur somme root2+root3 = -b-root1.
	
	complex<double> HSum= (-b/a-root1)/2.;				//demi somme des deux racines restantes.
	complex<double> Hdiff = sqrt(HSum*HSum + d/(a*root1));	//demi différence des deux racines restantes.
	//remarque : on a une division par root1. mais root1 peut être nul, par exemple pour le pol x^3. Dans ce cas, d est nul aussi...
	//le plus simple est d'intercepter ce comportement au début en testant d.
	
	roots[0]=root1;
	roots[1]=HSum+Hdiff;
	roots[2]=HSum-Hdiff;
	
	return 0;
};
*/