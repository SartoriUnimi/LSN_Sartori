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
#include "random.h"

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

vector<point> cities_circle;
vector<point> cities_square;
static const int Ncities = 34;

class Path{

    public:
        Path(){};
        ~Path(){};

        int M_Path[Ncities];
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
        void Sort(vector<point>&);
        int Select(Random &rnd);
        void Mutate(Random &rnd);
        int GetNpop() {return Npop;};
        double BestAv(vector<point>&);
        void crossing_over(Path);

    private:
        int Npop = 1000; //1000 ottimizza circle, 1800 ottimizza square
};

template <typename T> void Swap(T &a, T &b);


//Check function
int S_LEN = 0;
int N_LEN = 0;
int NGen = 500;
void Inizialisation(void);
int Pbc(int);
ofstream best_circle, av_circle, best_square, av_square, path_circle, path_square;



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