//métrique pour le cas P2L2.
//avec conditionnement

double Disc(double a, double b, double c, double d){
	return square(b)*square(c) - 4*a*cube(c) - 4*cube(b)*d + 18*a*b*c*d - 27*square(a)*square(d);
};
/*
//le calcul de la métrique. Probleme si le coef dominant est nul.
void MetricP2L2Unsafe(const double deriv[4], double metric[3]){
	double disc=Disc(deriv[0],3*deriv[1],3*deriv[2],deriv[3]);

	//cas dégénéré à gérer mieux...
	if(disc==0.){ metric[0]=1; metric[1]=0; metric[2]=1;
		return;
	};
	
	complex<double> phi[2][2];
	complex<double> roots[3];
	//si disc > 0, il y a une simplification, mais on pourra s'en occuper plus tard.
	
	//attention : dans le polynôme, les coefficients sont les dérivées multipliées par Cnk
	ThirdDegree(complex<double> (deriv[0],0.), complex<double> (3*deriv[1],0.), 
				complex<double> (3*deriv[2],0.), complex<double> (deriv[3],0.), 
				roots);
	
	double coef = deriv[0]*pow(2.*fabs(disc),-1/3.); coef *= disc<0?pow(2.,-1/6.):1; 
	const complex<double> deux = complex<double>(2.,0.);
	const complex<double> sqtrois = complex<double>(sqrt(3.),0.);
	const complex<double> isqtrois = complex<double>(0.,sqrt(3.));
	
	phi[0][0] = coef * ( roots[0]*(roots[1]+roots[2])-deux*roots[1]*roots[2] );
	phi[1][0] = coef * ( deux*roots[0]-(roots[1]+roots[2]) );
	
	phi[1][1] = coef *( roots[1]-roots[2] )*sqtrois;//*(disc>0?sqtrois:isqtrois);
	
	phi[0][1] = roots[0]*phi[1][1];
	
	//il faut inverser phi. On stocke toujours dans invphi, on utilise la formule tcom phi/det phi.
	//du moment que disc !=0, on n'est pas censé avoir de problème.
	complex<double> detphi = phi[0][0]*phi[1][1]-phi[0][1]*phi[1][0];
	complex<double> invphi[2][2];
	invphi[0][0] = phi[1][1]/detphi;
	invphi[1][1] = phi[0][0]/detphi;
	invphi[0][1] = -phi[0][1]/detphi;
	invphi[1][0] = -phi[1][0]/detphi;
	
	metric[0] = real(conj(invphi[0][0])*invphi[0][0] + conj(invphi[1][0])*invphi[1][0]);
	metric[1] = real(conj(invphi[0][0])*invphi[0][1] + conj(invphi[1][0])*invphi[1][1]);
	metric[2] = real(conj(invphi[0][1])*invphi[0][1] + conj(invphi[1][1])*invphi[1][1]);	
	
	//cette métrique définit une ellipse optimale dans {|pi|<=1}
	
	double detpow = pow(metric[0]*metric[2]-square(metric[1]), 1/8.);
	if(detpow!=0.){ //sinon on ne sait pas trop quoi faire.
	metric[0]/=detpow;
	metric[1]/=detpow;
	metric[2]/=detpow;
	};
	return;
};	
*/
inline double sign(double a){return a>0? 1: (a<0? -1 :0);};

double absPi(double a, double b, double c, double d, double z){ return fabs(a*cube(z) +b*square(z) +c*z+d)/cube(sqrt(square(z)+1));};

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
	//cas dégénéré à gérer mieux...
	if(disc==0.){ metric[0]=1; metric[1]=0; metric[2]=1;
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
	//if(CondSym2(metric)<=maxCond || maxCond<1){return;};
	
	//calcul du point le plus proche, par y Dx Pi - x Dy Pi
	double nearest[3]; //b, 2*c-a, d-2*b, -c
	if(ThirdDegreeReal(deriv[1], 2*deriv[2]-deriv[0], deriv[3]-2*deriv[1], -deriv[2], nearest) >0){// si cedans ce cas, trois points minima locaux, parmi lesquels if faut choisir.
		double dist[3];
		for(int i=0; i<3; i++){dist[i] = absPi(deriv[0], 3*deriv[1], 3*deriv[2], deriv[3], nearest[i]);};
		int pos = argmax(dist,3);
		nearest[0]=nearest[pos];
		//nearest[0] = nearest[1];
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
	
	if(maxCond==-2.){metric[0] = phi[0][0]; metric[1] = phi[0][1]; return;};
	if(maxCond==-3.){metric[0] = phi[1][0]; metric[1] = phi[1][1]; return;};
	if(maxCond==-4.){metric[0] = zPi[0];    metric[1] = zPi[1];    return;};
	if(maxCond==-5.){metric[0] = alpha;     metric[1] = beta;   metric[2] = lambdaC; return;};
	
	
	//zIPi[0] = invPhi[0][0]*zTPi[0] + invPhi[1][0]*zTPi[1]; zIPi[1] = invPhi[0][1]*zTPi[0] + invPhi[1][1]*zTPi[1];

	//ProdSquare2(phi,zIPi,zTPi);
	//metric[0]=zIPi[0]; metric[1]=zIPi[1];
	//return;
	
	/*metric[0]=alpha;
	metric[1]=beta;
	metric[2]=lambdaC;
	return;
	*/
	
	//conditionnement limite?
	double metricC[3];
	hLambda(invPhi, lambdaC, sdisc, metricC);
	double condC=CondSym2(metricC);
	
	//hLambda(invPhi, maxCond, sdisc, metricC);
	
	//metric[0]=metricC[0]; metric[1]=metricC[1]; metric[2]=metricC[2];
	//return;
	
	
	//hLambda(invPhi, 1./2., sdisc, metricC);	
	//metric[0]=metricC[0]; metric[1]=metricC[1]; metric[2]=metricC[2];
	//return;
	
	if(maxCond>0 && maxCond<condC){//ellispe liée à zPi. valeurs propres : lambdam et lambdaM
		double lambdam = pow(absPi(deriv[0], 3*deriv[1], 3*deriv[2], deriv[3], nearest[0]),2/3.);
		double lambdaM = lambdam/max(maxCond,1.);
				
		metric[0] = lambdam*square(zPi[0])+lambdaM*square(zPi[1]);
		metric[1] = (lambdam-lambdaM)*zPi[0]*zPi[1];
		metric[2] = lambdam*square(zPi[1])+lambdaM*square(zPi[0]);
		
	} else {//ellipse basée sur hLambda. 
		//metric[0]=metricC[0]; metric[1]=metricC[1]; metric[2]=metricC[2];
		//return;
		
		
		SquareToSym2_AAT(invPhi, metric); //on renverse l'ordre...
		
		if(maxCond == -6.){return;};
		
		//alpha = metric[0]-sdisc*metric[2]/3.;
		//beta = 4/3.*metric[2];
		alpha = metric[2]-sdisc*metric[0]/3.;
		beta = 4/3.*metric[0];
		double psiC = square(sqrt(maxCond)+1/sqrt(maxCond))*(metric[0]*metric[2]-square(metric[1]));
		double lambdaList[2]; 	//calcul du paramètre lambda présent, à l'aide de la borne sur le conditionnement
		SecondDegreeReal(3*square(alpha) +sdisc*psiC, 2*(3*alpha*beta-2*psiC), 3*square(beta), lambdaList);
		
		if(maxCond==-7){metric[0] = alpha; metric[1] = beta; metric[2] = psiC; return;};
		
		//metric[0] = lambdaList[0]; metric[1] = lambdaList[1]; metric[2] = psiC; return;
		
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

	if(true || abs(deriv[0]) == maxC){MetricP2L2UCond(deriv, metric, maxCond);};
/*
	double derivA[4];
	double metricA[3];

	if(abs(deriv[3]) == maxC){
		for(int i=0; i<4; i++){derivA[i] = deriv[3-i];};
		MetricP2L2Unsafe(derivA, metricA);
		for(int i=0; i<3; i++){metric[i] = metricA[2-i];};
	};

	if(abs(derivR[0]) == maxC){
		MetricP2L2Unsafe(derivR, metricA);
		metric[0] = ( metricA[0]+2*metricA[1]+metricA[2])/2;
		metric[1] = ( metricA[0]             -metricA[2])/2;
		metric[2] = ( metricA[0]-2*metricA[1]+metricA[2])/2;
	};

	if(abs(derivR[3]) == maxC){
		for(int i=0; i<4; i++){derivA[i] = derivR[3-i];};		
		MetricP2L2Unsafe(derivA, metricA);
		metric[0] = ( metricA[0]+2*metricA[1]+metricA[2])/2;
		metric[1] = (-metricA[0]             +metricA[2])/2;
		metric[2] = ( metricA[0]-2*metricA[1]+metricA[2])/2;
	};
	*/
	
	/*
	//facteur multiplicatif pour l'équilibrage de l'aire.
	double detpow = pow(metric[0]*metric[2]-square(metric[1]), 1/8.);
	if(detpow!=0.){ //sinon on ne sait pas trop quoi faire.
		metric[0]/=detpow;
		metric[1]/=detpow;
		metric[2]/=detpow;
	};
	 */
	
	return;
	//Des problèmes surviennent lorsque le coefficient dominant est trop petit.
	//idée : le polynôme ne peut pas être simultanément divisible par x,y,x+y et x-y.
	//On fait donc les rotations ssociées et on rend celle qui donne le coefficient dominant maximal.
};


 
