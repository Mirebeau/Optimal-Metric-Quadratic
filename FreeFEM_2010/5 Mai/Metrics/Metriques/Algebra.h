/*
 *  Algebra.h
 *  
 *
 *  Created by Jean-Marie Mirebeau on 04/06/09.
 *  Copyright 2009 UPMC. All rights reserved.
 *
 */


//la partie équations algébriques
complex<double> FirstDegree(complex<double> a, complex<double> b);
double FirstDegreeReal(double a, double b);
int SecondDegree(complex<double> a, complex<double> b, complex<double> c, complex<double> roots[]);
int SecondDegreeReal(double a, double b, double c, double roots[]);
int ThirdDegree(complex<double> a, complex<double> b, complex<double> c, complex <double> d, complex<double> roots[]);
int ThirdDegreeReal(double a, double b, double c, double d, double roots[]);
int ThirdDegreeRealDegen(double a, double b, double c, double d, double roots[]);


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
	if(delta<0){roots[0] = -b/(2*a); roots[1] = sqrt(-delta); return-1;};
	
	roots[0] = (-b-sqrt(delta))/(2*a);
	roots[1] = (-b+sqrt(delta))/(2*a);
	return 1;
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
		roots[2]=complex<double>(0.,0.);
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
inline double sqrt3(double x){return x>=0?pow(x,1./3.):-pow(-x,1/3.);}; //racine cubique avec signe
inline double maxAbs(double x, double y){return fabs(x)>fabs(y)?x:y;}; //revoie celui de plus grande valeur absolue
	
int ThirdDegreeReal(double a, double b, double c, double d, double roots[]){
	
	if(a==0.){roots[0]=TGV; return SecondDegreeReal(b,c,d,roots+1);};
	
	double alpha = -square(b)+3*a*c;
	double beta = -2*cube(b) +9*a*b*c - 27*square(a)*d;
	double delta = 4*cube(alpha) +square(beta);
	
	if(delta<=0){ //rappel : delta * disc<0
		complex<double> gamma = pow(complex<double>(beta, sqrt(-delta))/2.,1/3.); //rappel : éviter gamma = 0
		if(abs(gamma)==0.){roots[0]= sqrt3(-d/a); roots[1]=roots[0]; roots[2]=roots[0]; return 0;// c'est une racine triple (donc réelle)
		};
		complex<double> gamma2 = alpha/gamma;
		const complex<double> j = complex<double>(-0.5, sqrt(3.)/2);
		const complex<double> j2= complex<double>(-0.5, -sqrt(3.)/2);
		roots[0] = real(-b-gamma2+gamma)/(3*a);
		roots[1] = real(-b  -j*gamma2 +j2*gamma)/(3*a);
		roots[2] = real(-b -j2*gamma2 +j*gamma)/(3*a);
		return 1; //disc positif, trois racines réelles.
	} else {
		double sdelta = sqrt(delta);
		double gamma = sqrt3(maxAbs(beta+sdelta, beta-sdelta)/2.);
		if(abs(gamma)==0.){roots[0]= sqrt3(-d/a); roots[1]=roots[0]; roots[2]=roots[0]; return 0;// c'est une racine triple (donc réelle)
		};
		double gamma2 = alpha/gamma;
		roots[0] = (-b-gamma2+gamma)/(3*a);
		//roots[1]=gamma;
		//roots[2]=gamma2;
		//return -1;
		roots[1] = (-b +0.5*gamma2-0.5*gamma)/(3*a); //partie réelle
		roots[2] = (gamma2+gamma)/(2*sqrt(3)*a); //partie imaginaire
		return -1; //disc négatif, 1 réelle et deux complexes conjuguées.
	};
	return 0; //ne peut pas se produire
};


/*pas terrible.
il faut faire une routine qui exploite véritablement que 2 racines sont égales.
Par exemple en résolvant une eqn du second degré après avoir regardé les coefs. Ou bien avec PGCD P,P'
*/

/*Méthode : on résout P' = 0. Si racine double alors P a une racine triple.
 Sinon on identifie quelle est la bonne en injectant dans P et on retrouve l'autre 
 car on connait la somme*/

double evalThirdDegreeCircle(double a, double b, double c, double d, double r){
	double norm=sqrt(1+r*r);
	return (a*r*r*r +b*r*r +c*r +d)/(norm*norm*norm);
};
	
	
int ThirdDegreeRealDegen(double a, double b, double c, double d, double roots[]){
	double derRoots[3];
	int derSign = SecondDegreeReal(3*a, 2*b, c, derRoots);
	if(derSign<=0){roots[0]=derRoots[0]; roots[1] = derRoots[0]; roots[2] = derRoots[0]; return 3;}; //triple root
	
	double DRoot = (fabs(evalThirdDegreeCircle(a, b, c, d, derRoots[0]))< fabs(evalThirdDegreeCircle(a, b, c, d, derRoots[1])))? derRoots[0] : derRoots[1];
	
	double SRoot = -b/a - 2* DRoot;
	roots[0] = DRoot; roots[1] = DRoot; roots[2] = SRoot;
	return 2;
};
/*
const double CTripleRootPol = 0.05;

int ThirdDegreeRealDegen(double a, double b, double c, double d, double roots[]){
	
	if(a==0.){roots[0]=TGV; return SecondDegreeReal(b,c,d,roots+1);};
	
	double alpha = -square(b)+3*a*c;
	double beta = -2*cube(b) +9*a*b*c - 27*square(a)*d;
	double delta = 4*cube(alpha) +square(beta);
	
	if(delta<=0){ //rappel : delta * disc<0
		complex<double> gamma = pow(complex<double>(beta, sqrt(-delta))/2.,1/3.); //rappel : éviter gamma = 0
		if(abs(gamma)<=CTripleRootPol){roots[0]= sqrt3(-d/a); roots[1]=roots[0]; roots[2]=roots[0]; return 0;// c'est une racine triple (donc réelle)
		};
		complex<double> gamma2 = alpha/gamma;
		const complex<double> j = complex<double>(-0.5, sqrt(3.)/2);
		const complex<double> j2= complex<double>(-0.5, -sqrt(3.)/2);
		roots[0] = real(-b-gamma2+gamma)/(3*a);
		roots[1] = real(-b  -j*gamma2 +j2*gamma)/(3*a);
		roots[2] = real(-b -j2*gamma2 +j*gamma)/(3*a);
		return 1; //disc positif, trois racines réelles.
	} else {
		double sdelta = sqrt(delta);
		double gamma = sqrt3(maxAbs(beta+sdelta, beta-sdelta)/2.);
		if(abs(gamma)<=CTripleRootPol){roots[0]= sqrt3(-d/a); roots[1]=roots[0]; roots[2]=roots[0]; return 0;// c'est une racine triple (donc réelle)
		};
		double gamma2 = alpha/gamma;
		roots[0] = (-b-gamma2+gamma)/(3*a);
		//roots[1]=gamma;
		//roots[2]=gamma2;
		//return -1;
		roots[1] = (-b +0.5*gamma2-0.5*gamma)/(3*a); //partie réelle
		roots[2] = (gamma2+gamma)/(2*sqrt(3)*a); //partie imaginaire
		return -1; //disc négatif, 1 réelle et deux complexes conjuguées.
	};
	return 0; //ne peut pas se produire
};
*/
 
