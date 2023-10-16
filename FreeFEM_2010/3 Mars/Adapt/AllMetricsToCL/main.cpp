
#include  <iostream>
#include  <cfloat>
using namespace std;
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <complex>

/**** Des routines nécessaires ****/
#include "../Metriques/BasicMath.h"
#include "../Metriques/Mat22.h"
#include "../Metriques/Algebra.h"

/**** Les métriques ****/
#include "../Metriques/P2L2Degen.h"
#include "../Metriques/P2L2Cond.h"
#include "../Metriques/P1L2H1.h"

enum MType { P1L2, P1H1, P2L2, P2H1, P2L2degen, NotIdentified}; //type of the metric to be used.


/*format d'entrée : type de métrique,
 conditionnement maximal (ignoré si <1),
 dérivées
*/

int main (int argc, char * const argv[]) {
	
	if(argc<3){cout << "Error : not enough arguments. Usage : MetricType, Bound on conditionning, derivatives\n"; return 1;};
	
	char* CMetricType= argv[1];
	double maxCond = atof(argv[2]);
	double deriv[4];
	MType MetricType;
	
	if(!strcmp(CMetricType,"P1L2")){MetricType = P1L2;}
	else if(!strcmp(CMetricType,"P1H1")){MetricType = P1H1;}
	else if(!strcmp(CMetricType,"P2L2")){MetricType = P2L2;}
	else if(!strcmp(CMetricType,"P2H1")){MetricType = P2H1;}
	else if(!strcmp(CMetricType,"P2L2degen")){MetricType = P2L2degen;}
	else {MetricType = NotIdentified; cout << "Error : unidentified metric type\n"; return 1;}; 
	
	if( (MetricType == P1L2 || MetricType == P1H1) && argc!=6){
		cout << "Error, should be 5 arguments : MetricType, Bound on conditionning, order 2 derivatives."; return 1;
	};
	
	if(MetricType == P2L2 || MetricType == P2H1){
		if(argc!=7){cout << "Error, should be 6 arguments : MetricType, Bound on conditionning, order 3 derivatives."; return 1;};
		deriv[3] = atof(argv[6]);
	};
	
	deriv[0] = atof(argv[3]);
	deriv[1] = atof(argv[4]);
	deriv[2] = atof(argv[5]);
	deriv[3] = atof(argv[6]);
	
	
	double metric[3];
	if(MetricType==P1L2){MetricP1L2(deriv,metric, maxCond);};
	if(MetricType==P1H1){MetricP1H1(deriv,metric, maxCond);};
	if(MetricType==P2L2){MetricP2L2(deriv,metric, maxCond);};
	if(MetricType==P2H1){MetricP2H1(deriv,metric, maxCond);};
	if(MetricType==P2L2degen){MetricP2L2Degen(deriv, metric, maxCond);};
	
	ofstream metricfile;
	metricfile.open("MetricOut.txt");
	metricfile << metric[0] << "\n" << metric[1] << "\n" << metric[2] << "\n";
	
	//on renvoie qd même vers la command line
	
	cout << metric[0] << "\n" << metric[1] << "\n" << metric[2] << "\n";
	
	return 0;
}
