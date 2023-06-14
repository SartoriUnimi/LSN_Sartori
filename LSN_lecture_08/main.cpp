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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "SA.h"

using namespace std;

int main(){

  ofstream E;
  ofstream WriteConf;
  E.open("energy.dat");
  WriteConf.open("config.dat");

  //Inizializzazione
  Input();
  beta = 1./temp;

  //Annealing
  for(int i=1; i <= nSA; i++){
    cout << "SA block number " << i << endl;
    SAReset();                                //cancella i valori del passo precedente
    for(int iSA=1; iSA <= SAstep; iSA++){
      SAMove(nblk);                           //esegue la mossa SA
    }
    cout << "Acceptance rate " << SAaccepted/SAattempted << endl << endl;
    cout << "----------------------------" << endl << endl;
    E << beta << setw(12) << energy << setw(12) << err_energy << endl; //Salva le stime di energia per ogni passo
    WriteConf << i << setw(12) << mu << setw(12) << sigma << endl;  //Salva i valori di sigma e mu esplorati
    beta ++;
  }

  E.close();
  WriteConf.close();

  EnergyLevel_print(nblk);      //salva le stime di E per i migliori valori di sigma e mu

  ofstream sample;              //effettua un campionamento della funzione d'onda
  sample.open("psi.dat");
  for(int i = 0; i < 100000; i++){
    Move();
    sample << xold << endl;
  }

  sample.close();

  
  return 0;
}



void SAMove(int nblk){
  double p, m_old, s_old, energy_old, energy_new, err_old;

    // Mossa dell'annealing
  {

    //Old
      EnergyLevel(nblk);
      energy_old=energy;
      err_old=err_energy;
      m_old = mu;
      s_old = sigma;

    //New
      mu = m_old + deltapar*(rnd.Rannyu() - 0.5) ;
      sigma = s_old + deltapar*(rnd.Rannyu() - 0.5) ;
      EnergyLevel(nblk);

      energy_new = energy;

    //Metropolis test
      p =exp((energy_old-energy_new)*beta);

      if(p >= rnd.Rannyu())  
      {
      //Update
        m_old = mu;
        s_old = sigma;
        SAaccepted = SAaccepted + 1.0;
      } else {
        mu = m_old;
        sigma = s_old;
        energy = energy_old;
        err_energy = err_old;
      }
      SAattempted = SAattempted + 1.0;
      return;
  }
}

void EnergyLevel(int nblk){//Metropolis computation of Energy level
  for(int i = 0; i < 1000; i++) //equilibrate
    Move();
  for(int iblk=1; iblk <= nblk; iblk++) 
    {
    Reset(iblk);   //Reset block averages
    xold=0.;
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
  Averages(iblk);   //Save energy 
  }
}

void EnergyLevel_print(int nblk){//Metropolis computation of Energy level and prints
  for(int iblk=1; iblk <= nblk; iblk++) 
    {
    Reset(iblk);   //Reset block averages
    xold=0.;
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
  Averages_print(iblk);   //Save energy 
  }
}

void Input(void)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << "Double potential well (1-D)        " << endl;
  cout << "MC(NVT) simulation                 " << endl << endl;
  cout << "Potential v(x) = x^4 -5/2x^2" << endl << endl;
  cout << "Trial wave function Psi(x)=N*(exp(-(x-mu)^2/(2s^2))+exp(-(x+mu)^2/(2s^2)))" << endl << endl;
  cout << "The program uses natural units: h = 2*pi, m =1 " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input.in");

  Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> mu;
  cout << "Starting Mu = " << mu << endl;

  ReadInput >> sigma;
  cout << "Starting Sigma = " << sigma << endl;

  ReadInput >> delta;

  ReadInput >> temp;
  cout << "Starting Temperature = " << temp << endl;

  ReadInput >> deltapar;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> nSA;

  ReadInput >> SAstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  cout << "Number of SA blocks = " << nSA << endl;
  cout << "Number of steps in one block = " << SAstep << endl << endl;
  ReadInput.close();
  }

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
           glob_av = 0;
           glob_av2 = 0;
   }


   blk_av = 0;
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void SAReset (void) //Reset SA
{
  SAattempted = 0;
  SAaccepted = 0;
}

double psi2(double x,double m, double s){     //modulo quadro funzione d'onda

    return pow(exp(-(x-m)*(x-m)/(2.*s*s))+exp(-(x+m)*(x+m)/(2.*s*s)), 2);
}

double d2psi2(double x, double m, double s){  //parte cinetica
    return ((x-m)*(x-m)-s*s)/pow(s,4)*exp(-(x-m)*(x-m)/(2.*s*s))+((x+m)*(x+m)-s*s)/pow(s,4)*exp(-(x+m)*(x+m)/(2.*s*s));
}

void Move()
{
  double p, psi_old, psi_new;

    // Monte Carlo (NVT) move

    //Old
      psi_old = psi2(xold, mu, sigma);

    //New
      x = xold + delta*(rnd.Rannyu() - 0.5) ;

      psi_new = psi2(x, mu, sigma);

    //Metropolis test
      p =psi_new/psi_old;

      if(p >= rnd.Rannyu())  
      {
      //Update
        xold = x;
        accepted = accepted + 1.0;
      } else {
        x = xold;
      }
      attempted = attempted + 1.0;
      return;
}

void Measure() //Properties measurement
{

  walker = (x*x*x*x-2.5*x*x-0.5*d2psi2(x, mu, sigma)/sqrt(psi2(x, mu, sigma))); 

  return;
}

void Accumulate(void) //Update block averages
{
   blk_av = blk_av + walker;
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
    
    //cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;

    stima_etot = blk_av/blk_norm;
    glob_av += stima_etot;
    glob_av2 += stima_etot*stima_etot;
    err_etot=Error(glob_av,glob_av2,iblk);

//Salvataggio miglior stima
  if(iblk==nblk){
    energy = glob_av/(double)iblk;
    err_energy = err_etot;
    }
    //cout << "----------------------------" << endl << endl;

}

void Averages_print(int iblk) //Print results for current block
{
    ofstream Efin;
    Efin.open("energy_finale.dat",ios::app);
    //cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;

    stima_etot = blk_av/blk_norm;
    glob_av += stima_etot;
    glob_av2 += stima_etot*stima_etot;
    err_etot=Error(glob_av,glob_av2,iblk);

//Salvataggio miglior stima
  if(iblk==nblk){
    energy = glob_av/(double)iblk;
    err_energy = err_etot;
    }
  Efin << setw(12) << iblk <<  setw(12) << stima_etot << setw(12) << glob_av/(double)iblk << setw(12) << err_etot << endl;
    //cout << "----------------------------" << endl << endl;

  Efin.close();

}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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