//le calcul de la métrique P2H1.
 
void MetricP2H1(const double deriv[4], double metric[3]){
	double q1[3]; double q2[3]; //Formes quadratiques associées à Dx f et Dy f.
	q1[0] = deriv[0]; q1[1] = deriv[1]; q1[2] = deriv[2]; 
	q2[0] = deriv[1]; q2[1] = deriv[2]; q2[2] = deriv[3];
	
	double squareQ1[3]; double squareQ2[3];
	SquareSym2(q1,squareQ1); SquareSym2(q2,squareQ2);
	
	//On met la somme des carrés des matrices dans squareQ1
	squareQ1[0] += squareQ2[0]; squareQ1[1] += squareQ2[1]; squareQ1[2] += squareQ2[2];
	
	//la metrique est homogène à la racine carrée
	SqrtSym2(squareQ1,metric);
	
	double detpow = pow(metric[0]*metric[2]-square(metric[1]), 1/6.); //le det ne peut pas être négatif en théorie.
	if(detpow>0){
		metric[0]/=detpow; metric[1]/=detpow; metric[2]/=detpow;
	};	
	return;
};
