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
  dell'integrale di esercizio 1.1                            */

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
   ofstream g("data1.dat");

   // Definizione dei parametri del problema
   int M=100000;  // Numero totale di lanci
   int N=100;     // Numero di blocchi
   int L=M/N;     // Numero di lanci in ciascun blocco

   // Definizione dei vettori utilizzati per il calcolo
   vector<int> x(N);          // Indici dei blocchi
   vector<double> ave(N);     // Medie di ciascun blocco
   vector<double> av2(N);     // Quadrati delle medie di ciascun blocco
   vector<double> sum_prog(N); // Medie cumulate dei blocchi
   vector<double> su2_prog(N); // Quadrati delle medie cumulate dei blocchi
   vector<double> err_prog(N); // Errori sulle medie cumulate dei blocchi

   for(int i=1; i <= N; i++)
      x[i-1]=i;     //inizializzazione di x

   // Ciclo principale per il calcolo delle medie dei blocchi
   double sum=0.;
   for(int i=0; i<N; i++){
      sum = 0.;
      for(int j = 0; j < L; j++){
         double r = rnd.Rannyu(); // Generazione del numero casuale
         sum += r;
      }
      ave[i] = sum/L;    // Calcolo della media del blocco i-esimo
      av2[i] = ave[i]*ave[i]; // Calcolo del quadrato della media del blocco i-esimo
   }

   x=x*L;  // Restituisce il numero di lanci

   // Ciclo per il calcolo delle medie cumulate dei blocchi e degli errori
   for(int i=0; i<N; i++){
      for(int j=0; j<i+1; j++){
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i]=sum_prog[i]/(i+1);
      su2_prog[i]=su2_prog[i]/(i+1);
      err_prog[i]=error(sum_prog[i], su2_prog[i], i);

      // Scrittura dei dati su file di output
      g << setw(10) << x[i] << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;
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
