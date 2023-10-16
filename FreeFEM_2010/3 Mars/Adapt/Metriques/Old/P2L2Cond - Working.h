//métrique pour le cas P2L2.
//avec conditionnement

double Disc(double a, double b, double c, double d){
	return square(b)*square(c) - 4*a*cube(c) - 4*cube(b)*d + 18*a*b*c*d - 27*square(a)*square(d);
};

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

//le calcul de la métrique. Probleme si le coef dominant est nul.
void MetricP2L2UCond(const double deriv[4], double metric[3], double maxCond=-1.){
	double disc=Disc(deriv[0],3*deriv[1],3*deriv[2],deriv[3]);
	
	//cas dégénéré à gérer mieux...
	if(disc==0.){ metric[0]=1; metric[1]=0; metric[2]=1;
		return;
	};
	
	double phi[2][2];
	double roots[3];
	
	ThirdDegreeReal(deriv[0], 3*deriv[1], 3*deriv[2], deriv[3], roots);
		
	if(disc>0){// si trois racines réelles.
		phi[0][0] = roots[0]*(roots[1]+roots[2])-2.*roots[1]*roots[2];
		phi[1][0] = 2.*roots[0]-(roots[1]+roots[2]);
		phi[1][1] = sqrt(3.)*(roots[1]-roots[2]);
		phi[0][1] = roots[0]*phi[1][1];
	} else { // rappel : première racine réelle. puis partie réelle et imag des deux autres.
		phi[0][0] = roots[0]*2.*roots[1]-2.*(square(roots[1])+square(roots[2]));
		phi[1][0] = 2.*roots[0]-2.*roots[1];
		phi[1][1] = sqrt(3.)*2.*roots[2];
		phi[0][1] = roots[0]*phi[1][1];
	};
	
	double coef = deriv[0]*pow(2.*fabs(disc),-1/3.); coef *= disc<0?pow(2.,-1/6.):1; 
	phi[0][0]*=coef; phi[1][0]*=coef; phi[0][1]*=coef; phi[1][1]*=coef;
	
	double invPhi[2][2]; InvSquare2(phi,invPhi);
	SquareToSym2(invPhi,metric);
	
	double detpow = pow(metric[0]*metric[2]-square(metric[1]), 1/8.);
	if(detpow!=0.){ //sinon on ne sait pas trop quoi faire.
		metric[0]/=detpow;
		metric[1]/=detpow;
		metric[2]/=detpow;
	};
	return;
};	

void MetricP2L2(const double deriv[4], double metric[3]){
	//après rotation de pi/4. Attention : coef poly (degré 3)  = coef deriv *3 au milieu.
	double derivR[4];
	derivR[0] = (deriv[0]+3*deriv[1]+3*deriv[2]+deriv[3])/(2*sqrt(2));
	derivR[1] = (deriv[0]+  deriv[1]-  deriv[2]-deriv[3])/(2*sqrt(2));
	derivR[2] = (deriv[0]-  deriv[1]-  deriv[2]+deriv[3])/(2*sqrt(2));
	derivR[3] = (deriv[0]-3*deriv[1]+3*deriv[2]-deriv[3])/(2*sqrt(2));
	
	double maxC = max( max(abs(deriv[0]), abs(deriv[3])), max(abs(derivR[0]), abs(derivR[3])) );
	double derivA[4];
	double metricA[3];


	if(abs(deriv[0]) == maxC){MetricP2L2Unsafe(deriv, metric);return;};

	if(abs(deriv[3]) == maxC){
		for(int i=0; i<4; i++){derivA[i] = deriv[3-i];};
		MetricP2L2Unsafe(derivA, metricA);
		for(int i=0; i<3; i++){metric[i] = metricA[2-i];};
		return;
	};

	if(abs(derivR[0]) == maxC){
		MetricP2L2Unsafe(derivR, metricA);
		metric[0] = ( metricA[0]+2*metricA[1]+metricA[2])/2;
		metric[1] = ( metricA[0]             -metricA[2])/2;
		metric[2] = ( metricA[0]-2*metricA[1]+metricA[2])/2;
		return;
	};

	if(abs(derivR[3]) == maxC){
		for(int i=0; i<4; i++){derivA[i] = derivR[3-i];};		
		MetricP2L2Unsafe(derivA, metricA);
		metric[0] = ( metricA[0]+2*metricA[1]+metricA[2])/2;
		metric[1] = (-metricA[0]             +metricA[2])/2;
		metric[2] = ( metricA[0]-2*metricA[1]+metricA[2])/2;
		return;
	};
	
	// en cas d'erreur de max (?!?)
	metric[0] = 1; metric[1]=0; metric[2]=1;
	return;
	//Des problèmes surviennent lorsque le coefficient dominant est trop petit.
	//idée : le polynôme ne peut pas être simultanément divisible par x,y,x+y et x-y.
	//On fait donc les rotations ssociées et on rend celle qui donne le coefficient dominant maximal.
};


 
