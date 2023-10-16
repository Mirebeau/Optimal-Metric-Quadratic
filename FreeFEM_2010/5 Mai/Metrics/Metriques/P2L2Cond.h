//métrique pour le cas P2L2.
//avec conditionnement
const double CDoubleRoot = 1e-6;

double Disc(double a, double b, double c, double d){
	return square(b)*square(c) - 4*a*cube(c) - 4*cube(b)*d + 18*a*b*c*d - 27*square(a)*square(d);
};


double absPi(double a, double b, double c, double d, double z){ return fabs(a*cube(z) +b*square(z) +c*z+d)/cube(sqrt(square(z)+1));};

double distP1(double u, double v){
	double duv = fabs(atan(u) - atan(v));
	return min(duv, M_PI-duv); 
};

void orderP1(const double points[3], double ordered[3]){
	double dist[3]; 
	dist[0] = distP1(points[0],points[1]);
	dist[1] = distP1(points[1],points[2]);
	dist[2] = distP1(points[2],points[0]);
	
	int pos = argmin(dist, 3); //les deux plus proches doivent être en positions 1 et 2
	if(pos == 0)      {ordered[0] = points[2]; ordered[1] = points[1]; ordered[2] = points[0];}
	else if(pos == 1) {ordered[0] = points[0]; ordered[1] = points[1]; ordered[2] = points[2];}
	else              {ordered[0] = points[1]; ordered[1] = points[0]; ordered[2] = points[2];};//pos == 2
	return;
};

int hLambda(const double A[2][2], double lambda, double eps, double metric[3]){//invPhiPi^T Dlambda invPhiPi, A=invPhiPi
	//if(lambda<0 || lambda>(eps>0? 1: 2)){metric[0] = 1; metric[1] = 0; metric[2]=1; return 1;};
	double mu=(4-eps*cube(lambda))/(3*square(lambda));
	metric[0] = mu*square(A[0][0])+lambda*square(A[1][0]);
	metric[1] = mu*A[0][0]*A[0][1]+lambda*A[1][0]*A[1][1];
	metric[2] = mu*square(A[0][1])+lambda*square(A[1][1]);
	return 0;	
};

//le calcul de la métrique. Probleme si le coef dominant est nul.
void MetricP2L2UCond(const double deriv[4], double metric[3], double maxCond=-1.){
	double disc=Disc(deriv[0],3*deriv[1],3*deriv[2],deriv[3]);
	double sdisc = sign(disc);
	//cas dégénéré
	if(deriv[0]==0.){ metric[0]=1; metric[1]=0; metric[2]=1; return;};
	if(maxCond>=1. && fabs(disc)<=CDoubleRoot*square(square(deriv[0]))){ //square(square(deriv[0])+square(deriv[1])+square(deriv[2])+square(deriv[3]))
		metric[0]=1; metric[1]=0; metric[2]=1;
		return;
	};
	
	double phi[2][2];
	double roots[3];
	
	ThirdDegreeReal(deriv[0], 3*deriv[1], 3*deriv[2], deriv[3], roots);
		
	if(disc>0){// si trois racines réelles.
		//ordonnement des trois racines.
		double ordered[3];
		orderP1(roots, ordered);
		
		phi[0][0] = ordered[0]*(ordered[1]+ordered[2])-2.*ordered[1]*ordered[2];
		phi[1][0] = 2.*ordered[0]-(ordered[1]+ordered[2]);
		phi[1][1] = sqrt(3.)*(ordered[1]-ordered[2]);
		phi[0][1] = ordered[0]*phi[1][1];
	} else { // rappel : première racine réelle. puis partie réelle et imag des deux autres.
		phi[0][0] = roots[0]*2.*roots[1]-2.*(square(roots[1])+square(roots[2]));
		phi[1][0] = 2.*roots[0]-2.*roots[1];
		phi[1][1] = sqrt(3.)*2.*roots[2];
		phi[0][1] = roots[0]*phi[1][1];
	};
	
	double coef = deriv[0]*pow(2.*fabs(disc),-1/3.); //coef *= disc<0?pow(2.,-1/6.):1; 
	phi[0][0]*=coef; phi[1][0]*=coef; phi[0][1]*=coef; phi[1][1]*=coef;
	
	double invPhi[2][2]; InvSquare2(phi,invPhi);
	SquareToSym2(invPhi,metric);
	if(disc<0){metric[0]*=sqrt3(2.);metric[1]*=sqrt3(2.);metric[2]*=sqrt3(2.);};
	
	//calcul du conditionnement, s'il est bon, retour.
	if(CondSym2(metric)<=maxCond || maxCond<1){return;};
	
	//calcul du point le plus proche, par y Dx Pi - x Dy Pi
	double nearest[3]; //b, 2*c-a, d-2*b, -c
	if(ThirdDegreeReal(deriv[1], 2*deriv[2]-deriv[0], deriv[3]-2*deriv[1], -deriv[2], nearest) >0){// si cedans ce cas, trois points minima locaux, parmi lesquels if faut choisir.
		double dist[3];
		for(int i=0; i<3; i++){dist[i] = absPi(deriv[0], 3*deriv[1], 3*deriv[2], deriv[3], nearest[i]);};
		int pos = argmax(dist,3);
		nearest[0]=nearest[pos];
	};
	double zPi[2]; //destiné à contenir les deux composantes du point de ``fabs(pi)=1'' le plus proche de l'origine, mais normalisé.
	double normPi = sqrt(square(nearest[0])+1);
	zPi[0] = nearest[0]/normPi; zPi[1] = 1/normPi;
	
	//calcul du paramètre lambda critique.
	//identifié par le fait que zPi est vecteur propre.
	//double zPhiPi[2];    ProdSquare2(phi,    zPi,    zPhiPi);
	double zIPi[2]; ProdSquare2(invPhi, zPi, zIPi);
	double zTPi[2]; zTPi[0] = phi[0][0]*zPi[0] + phi[1][0]*zPi[1]; zTPi[1] = phi[0][1]*zPi[0] + phi[1][1]*zPi[1];
	
	//double alpha = zIPi[0]*zTPi[1] +sdisc*zIPi[1]*zTPi[0]/3.;
	//double beta  = 4*zIPi[1]*zTPi[0]/3.;
	double alpha = zIPi[1]*zTPi[0] +sdisc*zIPi[0]*zTPi[1]/3.;
	double beta  = 4*zIPi[0]*zTPi[1]/3.;
	double lambdaC; //lambda critique
	
	if(alpha == 0.){lambdaC = disc>0? 1 : sqrt3(2); //reprendre ce cas particulier...
	} else { lambdaC = sqrt3(beta/alpha);};
	
	//conditionnement limite?
	double metricC[3];
	hLambda(invPhi, lambdaC, sdisc, metricC);
	double condC=CondSym2(metricC);
		
	if(maxCond<condC){//ellispe liée à zPi. valeurs propres : lambdam et lambdaM
		double lambdam = pow(absPi(deriv[0], 3*deriv[1], 3*deriv[2], deriv[3], nearest[0]),2/3.);
		double lambdaM = lambdam/max(maxCond,1.);
		
		metric[0] = lambdam*square(zPi[0])+lambdaM*square(zPi[1]);
		metric[1] = (lambdam-lambdaM)*zPi[0]*zPi[1];
		metric[2] = lambdam*square(zPi[1])+lambdaM*square(zPi[0]);
	} else {
		SquareToSym2_AAT(invPhi, metric); //on renverse l'ordre...
		alpha = metric[2]-sdisc*metric[0]/3.;
		beta = 4/3.*metric[0];
		double psiC = square(sqrt(maxCond)+1/sqrt(maxCond))*(metric[0]*metric[2]-square(metric[1]));
		double lambdaList[2]; 	//calcul du paramètre lambda présent, à l'aide de la borne sur le conditionnement
		SecondDegreeReal(3*square(alpha) +sdisc*psiC, 2*(3*alpha*beta-2*psiC), 3*square(beta), lambdaList);
		lambdaList[1] = sqrt3(lambdaList[1]);
		if((lambdaC-lambdaList[1])*( (disc>0? 1: pow(2,1/3.)) - lambdaList[1])<0){lambdaList[0] = lambdaList[1];}
		else {lambdaList[0] = sqrt3(lambdaList[0]);};
		
		hLambda(invPhi, sqrt3(lambdaList[0]), sdisc, metric);
	};
	//retour de la matrice.
	
	return;
};	

void MetricP2L2(const double deriv[4], double metric[3], double maxCond=-1){
	//après rotation de pi/4. Attention : coef poly (degré 3)  = coef deriv *3 au milieu.
	double derivR[4];
	derivR[0] = (deriv[0]+3*deriv[1]+3*deriv[2]+deriv[3])/(2*sqrt(2));
	derivR[1] = (deriv[0]+  deriv[1]-  deriv[2]-deriv[3])/(2*sqrt(2));
	derivR[2] = (deriv[0]-  deriv[1]-  deriv[2]+deriv[3])/(2*sqrt(2));
	derivR[3] = (deriv[0]-3*deriv[1]+3*deriv[2]-deriv[3])/(2*sqrt(2));
	
	double maxC = max( max(abs(deriv[0]), abs(deriv[3])), max(abs(derivR[0]), abs(derivR[3])) );

	if(abs(deriv[0]) == maxC){MetricP2L2UCond(deriv, metric, maxCond);};

	double derivA[4];
	double metricA[3];

	if(abs(deriv[3]) == maxC){
		for(int i=0; i<4; i++){derivA[i] = deriv[3-i];};
		MetricP2L2UCond(derivA, metricA, maxCond);
		for(int i=0; i<3; i++){metric[i] = metricA[2-i];};
	};

	if(abs(derivR[0]) == maxC){
		MetricP2L2UCond(derivR, metricA, maxCond);
		metric[0] = ( metricA[0]+2*metricA[1]+metricA[2])/2;
		metric[1] = ( metricA[0]             -metricA[2])/2;
		metric[2] = ( metricA[0]-2*metricA[1]+metricA[2])/2;
	};

	if(abs(derivR[3]) == maxC){
		for(int i=0; i<4; i++){derivA[i] = derivR[3-i];};		
		MetricP2L2UCond(derivA, metricA, maxCond);
		metric[0] = ( metricA[0]+2*metricA[1]+metricA[2])/2;
		metric[1] = (-metricA[0]             +metricA[2])/2;
		metric[2] = ( metricA[0]-2*metricA[1]+metricA[2])/2;
	};
	
	
	
	//facteur multiplicatif pour l'équilibrage de l'aire.
	double detpow = pow(metric[0]*metric[2]-square(metric[1]), 1/8.);
	if(detpow!=0.){ //sinon on ne sait pas trop quoi faire.
		metric[0]/=detpow;
		metric[1]/=detpow;
		metric[2]/=detpow;
	};
	
	
	return;
	//Des problèmes surviennent lorsque le coefficient dominant est trop petit.
	//idée : le polynôme ne peut pas être simultanément divisible par x,y,x+y et x-y.
	//On fait donc les rotations ssociées et on rend celle qui donne le coefficient dominant maximal.
};


 
