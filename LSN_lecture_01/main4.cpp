/********************************************************************************
*********************************************************************************
    _/       _/_/_/_/  _/_/_/ Numerical Simulation Laboratory
   _/       _/       _/       Physics Department
  _/       _/  _/_/    _/     Universita' degli Studi di Milano
 _/       _/    _/       _/   Leonardo Giovanni Sartori
_/_/_/_/ _/_/_/_/  _/_/_/     email: leonardogiovanni.sartori@studenti.unimi.it
*********************************************************************************
*********************************************************************************/


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "VectorOperations.h"

using namespace std;

/*  Esercizio 1.2
    Produce i lanci dei 3 dadi e computa le medie su blocchi di lunghezza crescente*/
 
int main (int argc, char *argv[]){

   //inizializzazione generatore

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   //programma

   vector<string> tipi = {"standard", "esponenziale", "lorentziano"};
   vector<int> N = {1, 2, 10, 100}; 
   int n=10000;             

   for(int i = 0; i < tipi.size(); i++){
        for(int j = 0; j < N.size(); j++){
            ofstream g("Somme_"+to_string(N[j])+"punti_dado_"+tipi[i]+".dat");
            for(int k = 0; k < n; k++){
                double sum = 0.;
                for(int l = 0; l < N[j]; l++){
                    if(i==0)
                        sum += rnd.Rannyu();
                    if(i==1)
                        sum -= log(1-rnd.Rannyu());
                    if(i==2)
                        sum += tan(M_PI*(rnd.Rannyu()-0.5));
                }
                sum /= static_cast<double>(N[j]);
                g << sum << endl;
            }
            g.close();
        }
   }

   rnd.SaveSeed();
   return 0;
}

/********************************************************************************
*********************************************************************************
    _/       _/_/_/_/  _/_/_/ Numerical Simulation Laboratory
   _/       _/       _/       Physics Department
  _/       _/  _/_/    _/     Universita' degli Studi di Milano
 _/       _/    _/       _/   Leonardo Giovanni Sartori
_/_/_/_/ _/_/_/_/  _/_/_/     email: leonardogiovanni.sartori@studenti.unimi.it
*********************************************************************************
*********************************************************************************/
