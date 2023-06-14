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

/*Il seguente codice esegue il data blocking per la valutazione
  della devstd dell'integrale di esercizio 1.1               */

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

   ofstream g("data2.dat");
   int M=100000;  // Total number of throws
   int N=100;     // Number of blocks
   int L=M/N;     // Number of throws in each block
   vector<int> x(N);      
   for(int i=1; i <=N; i++)
      x[i-1]=i;     //[0,1,2,...,N-1]
   vector<double> ave(N);
   vector<double> av2(N);
   vector<double> sum_prog(N);
   vector<double> su2_prog(N);
   vector<double> err_prog(N);

   double sum=0.;
   for(int i=0; i<N; i++){
      sum = 0.;
      for(int j = 0; j < L; j++){
         double r = rnd.Rannyu(); // U[0,1) uniform distribution
         sum += (r-0.5)*(r-0.5);  // accumula la devstd
      }
      ave[i] = sum/L;
      av2[i] = ave[i]*ave[i];
   }

   x=x*L; //restituisce il numero di lanci

   //ciclo per la determinazione della media sui blocchi
   for(int i=0; i<N; i++){
      for(int j=0; j<i+1; j++){
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i]=sum_prog[i]/(i+1);
      su2_prog[i]=su2_prog[i]/(i+1);
      err_prog[i]=error(sum_prog[i], su2_prog[i], i);
      g << setw(10) << x[i] << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;
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

