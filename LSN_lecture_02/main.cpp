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
#include <string>
#include "random.h"
#include "IntegralMC.h"

using namespace std;
 
int main (int argc, char *argv[]){

   //controllo generatore

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
   ofstream g("integrale_unif.dat");
   ofstream g2("integrale_imp.dat");
   ofstream g3("integrale_bad.dat");
   int M = 100000;
   int N = 100;
   int L = M/N;

   double I_unif = 0.;  //accumulatore per l'integrale con campionamento uniforme
   double I2_unif = 0.;
   double I_imp = 0.;   //accumulatore per l'integrale con importance sampling
   double I2_imp = 0.;
   double I_bad = 0.;   //accumulatore per l'integrale con bad importance sampling
   double I2_bad = 0.;
   for(int i = 1; i <= N; i++){
      I_unif += IntegraUnif(rnd, L);
      I2_unif += Integra2Unif(rnd, L);
      double error = sqrt((I2_unif/i-I_unif/i*I_unif/i)/(i*L));
      g << i*L << "\t" << I_unif/i << "\t" << error << endl;

      I_imp += IntegraImportance(rnd, L);
      I2_imp += Integra2Importance(rnd, L);
      double errorimp = sqrt((I2_imp/i-I_imp/i*I_imp/i)/(i*L));
      g2 << i*L << "\t" << I_imp/i << "\t" << errorimp << endl;  

      I_bad += BadIntegraImportance(rnd, L);
      I2_bad += BadIntegra2Importance(rnd, L);
      double errorbad = sqrt((I2_bad/i-I_bad/i*I_bad/i)/(i*L));
      g3 << i*L << "\t" << I_bad/i << "\t" << errorbad << endl;     
   }


   g.close();
   g2.close();
   g3.close();
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
