/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <string>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

// far partire le simulazioni da alta t scendendo 
// e recuperare le configurazioni finali riduce la necessita' di equilibrare
int main()
{ 
  Input(); //Inizialization
  for(int i = temps; i >= 0; i--)
  {
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      /*for(int j = 0; j < 10000; j++)  //Equilibrate
        Move(metro);*/
      for(int istep=1; istep <= nstep; ++istep)
      {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk); //Compute averages of the current block
    }
    Print_values();   //Prints best values for each temperatures
    temp = temp - (2.0-0.5)/temps;
    beta = 1./temp;
  }
  // ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> temps;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int k;
  double cost;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
  k = (int)(rnd.Rannyu()*nspin);
	if (metro == 1){ // Metropolis sampling
      cost = 2 * s[k] * (J * (s[Pbc(k - 1)] + s[Pbc(k + 1)]) + h);
      attempted++;
      if (cost > 0) {
        if (rnd.Rannyu() < exp(-beta * cost)) {
          s[k] *= -1;
          accepted++;
        }
      } else {
        s[k] *= -1;
        accepted++;
      }
    } else { // Gibbs sampling
	attempted = 1;
	accepted = 1;
      s[k] = 1;
      cost = 2 * s[k] * (J * (s[Pbc(k - 1)] + s[Pbc(k + 1)]) + h);
      if (rnd.Rannyu() > (1. / (1. + std::exp(-beta * cost))))
        s[k] = -1;
    }
  }
}

/*double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}*/

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
    cout << "Block number " << iblk << " at Temperature " << temp << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;

    stima_c = beta*beta*(blk_av[ic]/blk_norm/(double)nspin-(double)nspin*stima_u*stima_u); //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;

    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;

    stima_x = blk_av[ix]/blk_norm*beta/(double)nspin; //Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;

    cout << "----------------------------" << endl << endl;
}

void Print_values(void)
{
  ofstream Ene, Heat, Mag, Chi;
  const int wd=12;

  string label;
  if(metro == 0)
    label = "-Gibbs";
  else
    label = "-Metro";
    
  //Energia interna
  Ene.open("output.ene"+label+".dat",ios::app);
  err_u=Error(glob_av[iu],glob_av2[iu],nblk);
  Ene << setw(wd) << temp <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
  Ene.close();

  //Capacità termica
  Heat.open("output.heat"+label+".dat",ios::app);
  err_c=Error(glob_av[ic],glob_av2[ic],nblk);
  Heat << setw(wd) << temp <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
  Heat.close();

  //Magnetizzazione
  Mag.open("output.mag"+label+".dat", ios::app);
  err_m=Error(glob_av[im],glob_av2[im],nblk);
  Mag << setw(wd) << temp <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
  Mag.close();

  //Suscettibilità
  Chi.open("output.chi"+label+".dat", ios::app);
  err_x=Error(glob_av[ix],glob_av2[ix],nblk);
  Chi << setw(wd) << temp <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
  Chi.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
