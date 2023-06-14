#include <cmath>
#include <vector>
#include "random.h"
using namespace std;

//Integra con sampling uniforme
double IntegraUnif(Random &rnd, int punti){ 
	double somma=0.;

	for(int i=0; i<punti; i++){
		double x= rnd.Rannyu();
		somma+=M_PI/2.*cos(M_PI*x/2.);
	}

	return somma/static_cast<double>(punti);
}

double Integra2Unif(Random &rnd, int punti){
	double somma=0.;

	for(int i=0; i<punti; i++){
		double x= rnd.Rannyu();
		somma+=pow(M_PI/2.*cos(M_PI*x/2.), 2);
	}

	return somma/static_cast<double>(punti);
}

//facciamo importance sampling con l'approssimazione di cos
double IntegraImportance(Random &rnd, int punti){
    double somma=0.;

	for(int i=0; i<punti; i++){
		double x= 1+sqrt(1-rnd.Rannyu());
		somma+=M_PI*cos(M_PI*x/2.)/(4.*(1-x));
	}

	return somma/static_cast<double>(punti);
}

double Integra2Importance(Random &rnd, int punti){
    double somma=0.;

	for(int i=0; i<punti; i++){
		double x= 1+sqrt(1-rnd.Rannyu());
		somma+=pow(M_PI*cos(M_PI*x/2.)/(4.*(1-x)), 2);
	}

	return somma/static_cast<double>(punti);
}

//facciamo importance sampling con sen
double BadIntegraImportance(Random &rnd, int punti){
    double somma=0.;

	for(int i=0; i<punti; i++){
		double x= 2./M_PI*acos(1-rnd.Rannyu());
		somma+=1./tan(M_PI*x/2.);
	}

	return somma/static_cast<double>(punti);
}
double BadIntegra2Importance(Random &rnd, int punti){
    double somma=0.;

	for(int i=0; i<punti; i++){
		double x= 2./M_PI*acos(1-rnd.Rannyu());
		somma+=pow(1./tan(M_PI*x/2.),2);
	}

	return somma/static_cast<double>(punti);
}