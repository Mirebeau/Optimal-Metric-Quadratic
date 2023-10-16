//les métriques pour les éléments P1, en norme L2 et H1.
//le calcul de la métrique P1L2.
//Ici il suffit de prendre la valeur absolue, multipliée par le déterminant à la bonne puissance (-1/3.).
int MetricP1L2(const double deriv[3], double metric[3], double maxCond=-1.){
	AbsSym2(deriv, metric);
	if(maxCond >= 0){
		double metricTemp[3];
		BoundExcentricity(metric, metricTemp, maxCond);
		metric[0] = metricTemp[0]; metric[1] = metricTemp[1]; metric[2] = metricTemp[2];
	};
	
	double det=metric[0]*metric[2]-square(metric[1]);
	if(det<=0){return 1;};
	double detpow = pow(det,-1/6.);
	metric[0] *= detpow;
	metric[1] *= detpow;
	metric[2] *= detpow;
	return 0;
};

//le calcul de la métrique P1H1.
//c'est le carré des dérivées, fois leur déterminant puissance (-1/2).
//le conditionnement maximal peut être précisé.
int MetricP1H1(const double deriv[3], double metric[3], double maxCond=-1.){
	double squareDeriv[3];
	SquareSym2(deriv, squareDeriv);
	if(maxCond>=0){BoundExcentricity(squareDeriv,metric,maxCond);
	} else { metric[0]=squareDeriv[0]; metric[1]=squareDeriv[1]; metric[2]=squareDeriv[2]; };
		
	double detpow=pow(metric[0]*metric[2]-square(metric[1]), 1/4.);
	if(detpow<=0){ return 0;};
	metric[0] /= detpow;
	metric[1] /= detpow;
	metric[2] /= detpow;
	return 0;
};


//le calcul de la métrique P2H1.

void MetricP2H1(const double deriv[4], double metric[3], double maxCond=-1.){
	double q1[3]; double q2[3]; //Formes quadratiques associées à Dx f et Dy f.
	q1[0] = deriv[0]; q1[1] = deriv[1]; q1[2] = deriv[2]; 
	q2[0] = deriv[1]; q2[1] = deriv[2]; q2[2] = deriv[3];
	
	double squareQ1[3]; double squareQ2[3];
	SquareSym2(q1,squareQ1); SquareSym2(q2,squareQ2);
	
	//On met la somme des carrés des matrices dans squareQ1
	squareQ1[0] += squareQ2[0]; squareQ1[1] += squareQ2[1]; squareQ1[2] += squareQ2[2];
	
	//la metrique est homogène à la racine carrée
	SqrtSym2(squareQ1,metric);

	if(maxCond >= 0){
		double metricTemp[3];
		BoundExcentricity(metric, metricTemp, maxCond);
		metric[0] = metricTemp[0]; metric[1] = metricTemp[1]; metric[2] = metricTemp[2];
	};
	
	double detpow = pow(metric[0]*metric[2]-square(metric[1]), 1/6.); //le det ne peut pas être négatif en théorie.
	if(detpow>0){
		metric[0]/=detpow; metric[1]/=detpow; metric[2]/=detpow;
	};	
	return;
};
