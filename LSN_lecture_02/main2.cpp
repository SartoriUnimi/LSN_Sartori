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
#include <vector>
#include <cmath>
#include "random.h"
#include "VectorOperations.h"

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
    ofstream g("rw_cubico.dat");
    ofstream h("rw_continuo.dat");
    int n = 10000; //number of walks
    int M = 200; //number of steps
    int N = 100; //number of blocks
    int L = n/N; //walks per block


    vector<double> sum2_cubic(M);       //accumulatori dei valori quadrati sugli M passi
    vector<double> sum2_continuum(M);
    vector<double> avg_cubic(M);        //accumulatori delle medie sugli M passi
    vector<double> avg_continuum(M);

    for(int k = 1; k <= N; k++){
        vector<double> sum_cubic(M);    //accumulatori delle somme sugli M passi
        vector<double> sum_continuum(M);
        for(int l = 0; l < L; l++){
            vector<double> rw_cubic(3);     //rappresentano il RW in coordinate 3D
            vector<double> rw_continuum(3);
            for(int j=0; j < M; j++){
                double r = rnd.Rannyu();
                //step di RW cubico
                int step = 0;
                if((int)floor(6*r)%2==0)            //sceglie verso
                    step ++;
                else
                    step --;
                rw_cubic[(int)floor(6*r)%3] += step;//e direzione

                //step di RW continuo
                double phi = 2*M_PI*r;              //sceglie azimuth
                double theta = M_PI * rnd.Rannyu(); //sceglie elevazione
                vector<double> stepc = {sin(phi)*cos(theta), sin(theta)*sin(phi), cos(phi)};
                rw_continuum +=  stepc; 

                sum_cubic[j]+=sqrt(rw_cubic*rw_cubic);
                sum_continuum[j]+=sqrt(rw_continuum*rw_continuum);
            }
        }
        avg_cubic += sum_cubic/(double)L;           //calcola le medie
        avg_continuum += sum_continuum/(double)L;

        for(int j = 0;  j < M; j++){                //i quadrati delle medie
        sum2_cubic[j]+= sum_cubic[j]*sum_cubic[j]/(L*L);
        sum2_continuum[j]+= sum_continuum[j]*sum_continuum[j]/(L*L);
        }
    }
    
    //scrive su file
    for(int i = 0; i < M; i++){
        g << (i+1) << "\t" << avg_cubic[i]/N << "\t" << sqrt((sum2_cubic[i]/N-pow(avg_cubic[i]/N, 2))/N) << endl;
        h << (i+1) << "\t" << avg_continuum[i]/N << "\t" << sqrt((sum2_continuum[i]/N-pow(avg_continuum[i]/N, 2))/N) << endl;
    }

   g.close();
   h.close();
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
