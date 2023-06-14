/********************************************************************************
*********************************************************************************
    _/       _/_/_/_/  _/_/_/ Numerical Simulation Laboratory
   _/       _/       _/       Physics Department
  _/       _/  _/_/    _/     Universita' degli Studi di Milano
 _/       _/    _/       _/   Leonardo Giovanni Sartori
_/_/_/_/ _/_/_/_/  _/_/_/     email: leonardogiovanni.sartori@studenti.unimi.it
*********************************************************************************
*********************************************************************************/

#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
double vtail, ptail, bin_size, nbins, sd;
double walker;

// averages
double blk_av, blk_norm, accepted, attempted, SAaccepted, SAattempted;
double glob_av, glob_av2;
double stima_etot;
double err_etot;

//configuration
double x;
double xold;


// simulation
int nstep, nblk, nSA, SAstep;
double delta, temp, beta, energy, err_energy, deltapar;

//psi parameters
double mu, sigma;

//functions
void Input(void);
void Reset(int);
void SAReset(void);
void Accumulate(void);
void Averages(int);
void Averages_print(int);
void Move(void);
void SAMove(int);
void Measure(void);
void EnergyLevel(int);
void EnergyLevel_print(int);
double psi2(double, double, double);
double d2psi2(double, double, double);
double Error(double,double,int);

#endif

/********************************************************************************
*********************************************************************************
    _/       _/_/_/_/  _/_/_/ Numerical Simulation Laboratory
   _/       _/       _/       Physics Department
  _/       _/  _/_/    _/     Universita' degli Studi di Milano
 _/       _/    _/       _/   Leonardo Giovanni Sartori
_/_/_/_/ _/_/_/_/  _/_/_/     email: leonardogiovanni.sartori@studenti.unimi.it
*********************************************************************************
*********************************************************************************/
