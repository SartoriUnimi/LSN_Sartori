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
#include <cmath>
#include "random.h"


using namespace std;
 
int main (int argc, char *argv[]){

   //controllo generatore

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
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
   }
      
   ofstream gC1("call_direct.dat");
   ofstream gP1("put_direct.dat");
   ofstream gC2("call_discretize.dat");
   ofstream gP2("put_discretize.dat");

   double S0 = 100.;    //asset price a t = 0
   double T = 1.;       //delivery time
   double K =100.;      //strike prize
   double r = 0.1;      //risk-free interest rate
   double sigma = 0.25; //volatility

   double M = 100000;    //number of tries
   double N = 100;      //number of blocks
   double L = M/N;      //tries per block
   double n = 100;      //number of intervals

   double sumC = 0.;    //accumulatori per il caso diretto
   double sumC2 = 0.;
   double sumP = 0.;
   double sumP2 = 0.;

   double sumc = 0.;    //accumulatori per il caso discreto
   double sumc2 = 0.;
   double sump = 0.;
   double sump2 = 0.;

   for(int i = 1; i <= N; i++){
    for(int j = 0; j < L; j++){

        //direct
        double W =  rnd.Gauss(0., 1.);
        double St = S0 * exp((r - sigma*sigma/2.)*T + sigma*W);
        double C = 0.;
        double P = 0.; 
        if (St - K > 0) C += exp(-r*T)*(St-K);
        else P -= exp(-r*T)*(St-K);
        sumC += C;
        sumC2 += C*C;
        sumP += P;
        sumP2 += P*P;

        //discretized
        double st = S0;
        for(int k = 1;  k <= n; k++){
            double z = rnd.Gauss(0., 1.);
            st = st * exp((r - sigma*sigma/2.)*T/n + sigma*z*sqrt(T/n));
        }
        double c = 0.;
        double p = 0.;
        if (st - K > 0) c += exp(-r*T)*(st-K);
        else p -= exp(-r*T)*(st-K);
        sumc += c;
        sumc2 += c*c;
        sump += p;
        sump2 += p*p;
    }

   // calcolo di medie e errori
    double avgC = sumC/(i*L);
    double avgC2 = sumC2/(i*L);
    double errorC = sqrt((avgC2 - avgC*avgC)/(i*L));
    double avgP = sumP/(i*L);
    double avgP2 = sumP2/(i*L);
    double errorP = sqrt((avgP2 - avgP*avgP)/(i*L));
    double avgc = sumc/(i*L);
    double avgc2 = sumc2/(i*L);
    double errorc = sqrt((avgc2 - avgc*avgc)/(i*L));
    double avgp = sump/(i*L);
    double avgp2 = sump2/(i*L);
    double errorp = sqrt((avgp2 - avgp*avgp)/(i*L));

    //scrittura su file
    gC1 << i << "\t" << avgC << "\t" << errorC << endl;
    gP1 << i << "\t" << avgP << "\t" << errorP << endl;
    gC2 << i << "\t" << avgc << "\t" << errorc << endl;
    gP2 << i << "\t" << avgp << "\t" << errorp << endl;
   }


   gC1.close();
   gP1.close();
   gC2.close();
   gP2.close();
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
