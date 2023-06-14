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
#include "random.h"
#include "Functions.h"
#include "VectorOperations.h"

using namespace std;

/* Effettua il test del chi2 sulla distribuzione generata 
   per il calcolo dell'integrale dell'esercizio 1.1      */
 
int main (int argc, char *argv[]){

   //inizializzazione generatore

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
      Primes >> p1 >> p2 ;
      Primes >> p1 >> p2 ;
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

   ofstream g("data3.dat");
   int M=100;               // # of bin, # of subintervals
   int n=1000000;             // # of throws in each block
   vector<double> chi2(M);
   vector<int> bin(M);      // occupazione di ciascun bin

   // calcola chi2
   for(int k=0; k<M; k++){
        for(int i=0; i<n; i++){
            double r = rnd.Rannyu();
            bin[floor(r*M)]++;
        }
        for(int i=0; i<M; i++){
            chi2[k] += pow(bin[i]-n/M,2)/(n/M);
            bin[i]=0;
        }
        g << chi2[k] << endl;
   }

   g.close();

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
