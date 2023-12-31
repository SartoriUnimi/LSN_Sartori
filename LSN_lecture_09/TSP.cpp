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
#include <iomanip>
#include "random.h"
#include "VectorOperations.h"
#include "TSP.h"

using namespace std;

int main(){

    //inizializzazione

    Inizialisation();

    Population pop_circle, pop_square;
    pop_circle.Birth(rnd);
    pop_square.Birth(rnd);
    pop_circle.Sort(cities_circle);
    pop_square.Sort(cities_square);
    for(int i = 0; i < NGen; i++)
    {
        cout << "===========================================\n" << endl;
        cout << "Generazione : " << i+1 << endl << endl;

        pop_circle.Mutate(rnd);
        pop_square.Mutate(rnd);
        pop_circle.Sort(cities_circle);
        pop_square.Sort(cities_square);

        best_circle << i+1 << setw(12) << pop_circle.Chr[0].L2(cities_circle)<< endl;
        best_square << i+1 << setw(12) << pop_square.Chr[0].L2(cities_square)<< endl;
        av_circle << i+1 << setw(12) << pop_circle.BestAv(cities_circle)<< endl;
        av_square << i+1 << setw(12) << pop_square.BestAv(cities_square)<< endl;
    }


    pop_circle.Chr[0].Print_cities(cities_circle, path_circle);
    pop_square.Chr[0].Print_cities(cities_square, path_square);

    best_circle.close();
    best_square.close();
    av_circle.close();
    av_square.close();
    path_circle.close();
    path_square.close();

    return 0;
}

void Inizialisation(void)
{
    best_circle.open("circle_best.dat");
    av_circle.open("circle_average.dat");
    best_square.open("square_best.dat");
    av_square.open("square_average.dat");
    path_circle.open("best_final_circle.dat");
    path_square.open("best_final_square.dat");

    //Read seed for random numbers
    ifstream Primes, Seed;
    int p1, p2;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    for(int i = 0; i < Ncities; i++)
    {   
        S_LEN += i+1;                             //per il check
        N_LEN += (i+1)*(i+1);
        double a = rnd.Rannyu()*2*M_PI;                 
        point q(cos(a), sin(a));
        cities_circle.push_back(q);             //cities on circumference, assegna la posizione angolare
        point p(rnd.Rannyu(), rnd.Rannyu());
        cities_square.push_back(p);             //cities on square, assegna posizione cartesiana
    }
}

bool Path :: Check()
{
    if(M_Path[0] != 1) {                                        //checking the first position
        //Print();
        return false;
        }   
    int Norm = 0;                                               //checking no cities echoed
    int Sum = 0;
    for(int i = 0; i < Ncities; i++)
    {
        Sum += M_Path[i];
        Norm += M_Path[i]*M_Path[i];
    }
    if(Sum != S_LEN) {                                          //with sum
        //Print();
        return false;
        }
    if(Norm != N_LEN) {                                         //with norm
        //Print();
        return false;
        }
    else
        return true;
}

void Path :: Empty()
{
    for(int i = 0; i < Ncities; i++)
        M_Path[i] = 0;
}

void Path :: Fill(Random &rnd)
{
    if(M_Path[0] != 0) {this->Empty();};
    int p = floor(rnd.Rannyu()*(Ncities-1)+1);
    M_Path[0] = 1;
    for(int i = 2; i<= Ncities; i++)
    {
        while(M_Path[p] != 0)
            p = floor(rnd.Rannyu()*(Ncities-1)+1);
        M_Path[p] = i;
    }
}


void Population :: Birth(Random &rnd)
{
    if(Chr.size() != 0) {Chr.clear();}
    for(int i = 0; i < Npop; i++)
    {
        Path p;
        p.Fill(rnd);
        Chr.push_back(p);
    }
}

double Path :: L2(vector<point>& cities)
{
    double len = 0.;
    for(int i = 0; i < Ncities-1; i++){
        len += pow(cities[M_Path[i]-1].x()-cities[M_Path[i+1]-1].x(),2) + pow(cities[M_Path[i]-1].y()-cities[M_Path[i+1]-1].y(),2);
    }
            len += pow(cities[M_Path[Ncities-1]-1].x()-cities[M_Path[0]-1].x(),2) + pow(cities[M_Path[Ncities-1]-1].y()-cities[M_Path[0]-1].y(),2);
    return len;
}

void Population :: Sort(vector<point>& cities)
{
    for(int i = 0; i < Npop - 1; i++){
        for(int j = i+1; j < Npop; j++){
            if(Chr[i].L2(cities) > Chr[j].L2(cities))
                Swap(Chr[i], Chr[j]);
        }
    }
    /*for(int i = 0; i < Npop; i++)
        cout  << Chr[i].L2(cities) << endl;
    cout << endl;*/
}

int Population :: Select(Random &rnd)
{
    return (int)((Npop)*pow(rnd.Rannyu(), 2.0)); //sistemare l'esponente
}

void Population :: Mutate(Random &rnd)
{
    Path Chr_appo[Npop];
    for(int i = 0; i < Npop; i++)
    {
        int ap = Select(rnd);
        Chr_appo[i] = Chr[ap];
    }
    for(int p = 0; p < Npop; p++)
    {
        double prob = rnd.Rannyu();
        if(prob < 0.4)                //pair permutation
            Chr_appo[p].pair_permutation();
        if(Chr_appo[p].Check()==false){throw runtime_error("Permutation Misfunction");}

        prob = rnd.Rannyu();
        if(prob < 0.2)                //shift of n positions of m continguous cities  
            Chr_appo[p].n_shift();
        if(Chr_appo[p].Check()==false){throw runtime_error("Shift Misfunction");}

        prob = rnd.Rannyu();
        if(prob < 0.2)                //permutation among m contiguous cities with other m different cities   
            Chr_appo[p].m_permutation();
        if(Chr_appo[p].Check()==false){throw runtime_error("PermutaM Misfunction");}

        prob = rnd.Rannyu();
        if(prob < 0.1)                //inversion of m cities
            Chr_appo[p].inversion();
        if(Chr_appo[p].Check()==false){throw runtime_error("Inversion Misfunction");}

        prob = rnd.Rannyu();
        if(prob < 0.7)                  //crossing over
            crossing_over(Chr_appo[p]);
    }
    for(int i = 0; i < Npop; i++)
    {
        Chr[i] = Chr_appo[i];
    }

}

double Population :: BestAv(vector<point>& cities){
    double sum = 0.;
    for(int i = 0; i < (int)(Npop/2); i++)
        sum += Chr[i].L2(cities);
    sum /= (int)(Npop/2);

    return sum;
}

void Path :: Print(){
    for(int i = 0; i < Ncities; i++)
        cout << M_Path[i] << "\t";
    cout << endl;
}

Path& Path :: operator= (const Path& pth){
    for(int i = 0; i < Ncities; i++)
        M_Path[i] = pth.M_Path[i];

    return *this;
}

void Path :: Print_cities(vector<point>& cities, ofstream &g){
    for(int i = 0; i < Ncities; i++)
        g << cities[M_Path[i]-1].x() << "\t" << cities[M_Path[i]-1].y() << "\n";
}

template <typename T> void Swap(T &a, T &b)
{
    T dep;
    dep = a;
    a = b;
    b = dep;
}

int Pbc(int pos)
{
    if(pos>(Ncities-1) && pos%(Ncities-1)!=0)    pos = pos%(Ncities-1); // non Ncities, devo saltare la posizione 0
    if(pos>(Ncities-1) && pos%(Ncities-1)==0)    pos = Ncities-1;
    if(pos < 0)       pos = Ncities + pos%Ncities;
    return pos;
}

void Path :: pair_permutation()
{
    //Print();
    //cout << "permuta coppie\n";
    int k = floor((Ncities-1)*rnd.Rannyu()+1);
    int l = floor((Ncities-1)*rnd.Rannyu()+1);
    Swap(M_Path[k], M_Path[l]);
    //Print();
}

void Path :: n_shift()
{
    //Print();
    //cout << "shift\n";
    int a = floor((Ncities-1)*rnd.Rannyu()+1);
    int n = floor((Ncities-1)*rnd.Rannyu()+1);
    int m = floor((Ncities-1)*rnd.Rannyu()+1);
    for(int i = 0; i < m; i++)
        Swap(M_Path[Pbc(a+i)], M_Path[Pbc(a+n+i)]);
    //Chr[p].Print();
}

void Path :: m_permutation()
{
    //cout << "PermutaM\n";
    //Print();
    int a = floor((Ncities-1)*rnd.Rannyu()+1);
    int b = floor((Ncities-1)*rnd.Rannyu()+1);
    int m = floor((Ncities-1)*rnd.Rannyu()/2.+1);
    for(int i = 0; i < m; i++)
        Swap(M_Path[Pbc(a+i)], M_Path[Pbc(b+i)]);
    //Print();
}

void Path :: inversion(){
    //cout << "inversione\n";
    //Print();
    int a = floor((Ncities-1)*rnd.Rannyu()+1);
    int m = floor((Ncities-1)*rnd.Rannyu()/2.);
    for(int i = 0; i < m; i++)
        Swap(M_Path[Pbc(a+i)], M_Path[Pbc(a+m-i)]);
    //Print();
}

void Population :: crossing_over(Path Chr_appo)
{
    //cout << "crossing over\n";
    int cut = floor((Ncities-1)*rnd.Rannyu()+1);
    int len = floor((Ncities-1)*rnd.Rannyu()+1);
    int a[len];
    int cont = 0;
    int m = Select(rnd);
    //Chr[m].Print();
    //Chr_appo.Print();
    for(int j = 0; j < len; j++)
    {
        a[j] = Chr_appo.M_Path[Pbc(cut+j)];
    } 

    for(int i  = 1; i < Ncities; i++)
    {
        for(int j = 0; j < len; j++)
        {
            if(Chr[m].M_Path[i]==a[j])
            {
                Chr_appo.M_Path[Pbc(cut+cont)]=Chr[m].M_Path[i];
                cont ++;
            }
        }
    }
    //Chr[m].Print();
    //Chr_appo.Print();
    if(Chr[m].Check()==false){throw runtime_error("Crossing-over Misfunction");}
    if(Chr_appo.Check()==false){throw runtime_error("Crossing-over Misfunction f");}
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