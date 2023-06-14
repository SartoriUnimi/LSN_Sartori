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

/* Esercizio 1.3
   Esegue la simulazione dell'esperimento di Buffon */
 
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

   // Apertura del file di output
   ofstream g("buffon.dat");

    // Definizione dei parametri
    const double L = 0.1;
    const double d = 0.2;
    const int n_trials = 1000000; // numero totale di lanci
    const int n_blocks = 100; // numero di blocchi
    const int n_steps = n_trials / n_blocks; // numero di lanci per blocco

   // Ciclo sui blocchi
    double sum = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < n_blocks; i++) {
        // Ciclo sui lanci nel blocco
        int count = 0;
        for (int j = 0; j < n_steps; j++) {
            count += lancio(L, d, rnd);
        }
        double pi_estimate = 2.0 * L * n_steps / (d * count);
        sum += pi_estimate;
        sum2 += pi_estimate * pi_estimate;

        // Scrittura dei risultati su file
        double err = error(sum/(i+1), sum2/(i+1), i);
        double avg = sum / (i + 1);
        g << (i + 1) * n_steps << "\t" << fixed << setprecision(6) << avg << "\t" << err << endl;
    }

   // Chiusura del file di output
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
