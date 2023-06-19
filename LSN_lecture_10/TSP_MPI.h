/********************************************************************************
*********************************************************************************
    _/       _/_/_/_/  _/_/_/ Numerical Simulation Laboratory
   _/       _/       _/       Physics Department
  _/       _/  _/_/    _/     Universita' degli Studi di Milano
 _/       _/    _/       _/   Leonardo Giovanni Sartori
_/_/_/_/ _/_/_/_/  _/_/_/     email: leonardogiovanni.sartori@studenti.unimi.it
*********************************************************************************
*********************************************************************************/

#ifndef __TSP__
#define __TSP__

#include <vector>
#include <fstream>
#include "random.h"
#include "mpi.h"
using namespace std;

//Random numbers
int seed[4];
Random rnd;

class point
{
public:

    point(){};
    ~point(){};

    point(double a,  double b)
    {
        m[0] = a;
        m[1] = b;
    }
    double x() {return m[0];};
    double y() {return m[1];};

private:    
    double m[2];
};

vector<point> cities;
static const int Ncities = 50;

class Path{

    public:
        Path(){};
        ~Path(){};

        int M_Path[Ncities];
        double beta = 0.;
        void Empty();
        void Fill(Random &rnd);
        bool Check();
        void Print();
        double L2(vector<point>&);
        Path& operator= (const Path&);
        void Print_cities(vector<point>&, ofstream &);

        void pair_permutation();
        void n_shift();
        void m_permutation();
        void inversion();
    
};

class Population{
    
    public:
        Population(){};
        ~Population(){};

        vector<Path> Chr;
        void Birth(Random &rnd);
        void Birth(Random &rnd, int);
        void Sort(vector<point>&);
        int Select(Random &rnd);
        void Mutate(Random &rnd, vector<point>&);
        int GetNpop() {return Npop;};
        double BestAv(vector<point>&);
        void crossing_over(Path);

    private:
        int Npop = 30; 
};

template <typename T> void Swap(T &a, T &b);


//Check function
int S_LEN = 0;
int N_LEN = 0;
int nSA = 5000;
int nblk = 1000;
int SAattempted, SAaccepted;
void Initialisation(ofstream&, int);
void SAReset(void);
int Pbc(int);
void Adjacent_Exchange(Path, Path);
void Adjacent_Exchange_MPI(int &rank, int &size, Population pop, MPI_Status stat);
ofstream av, path;



#endif

/********************************************************************************
*********************************************************************************
    _/       _/_/_/_/  _/_/_/ Numerical Simulation Laboratory
   _/       _/       _/       Physics Department
  _/       _/  _/_/    _/     Universita' degli Studi di Milano
 _/       _/    _/       _/   Leonardo Giovanni Sartori
_/_/_/_/ _/_/_/_/  _/_/_/     email: leonardogiovanni.sartori@studenti.unimi.it
*********************************************************************************
*********************************************************************************/
