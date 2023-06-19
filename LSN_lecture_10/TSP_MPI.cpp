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
#include "TSP_MPI.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]){

    // ====== INIZIALIZZAZIONE MPI ======
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat[size];
    // ====== INIZIALIZZAZIONE PROBLEMA ======
    ofstream best[size];
    Initialisation(best[rank], rank);

    //per ogni core costruiamo population di nblk elementi
    Population pop[size];
    pop[rank].Birth(rnd, nblk);
    pop[rank].Sort(cities);
    for(int i = nblk; i > 0; i--)
    {
	    pop[rank].Chr[nblk - i].beta = i*0.1*pow(0.001, log(rank + 2));
    }

    for(int iSA = 1; iSA <= nSA; iSA++)
    {
        SAReset();
        // Mutiamo le popolazioni
        pop[rank].Mutate(rnd, cities);
        MPI_Barrier(MPI_COMM_WORLD);

        //implementiamo lo scambio tra temperature adiacenti
        //nello stesso nodo
        for(int i = 0; i < nblk-1; i++)
            Adjacent_Exchange(pop[rank].Chr[i], pop[rank].Chr[i+1]);
        
        // tra nodi adiacenti
        Adjacent_Exchange_MPI(rank, size, pop[rank], stat[rank]);

        //aggiorniamo le temperature
        for(int jl = 0; jl < nblk; jl++)
        {
            pop[rank].Chr[jl].beta *= 1.000001;
        }
        MPI_Barrier(MPI_COMM_WORLD);


        pop[rank].Sort(cities);
        cout << "Rank " << rank << ": SA block number " << iSA << endl;
        cout << "Acceptance rate " << (double)(SAaccepted)/(double)(SAattempted) << endl << endl;
        cout << "----------------------------" << endl << endl;
        best[rank] << iSA << setw(12) << pop[rank].Chr[0].L2(cities)  << endl;
        if(rank == 0) //la migliore stima è alla temperatura più bassa, nodo 0, individuo 0
        {
            av << iSA << setw(12) << pop[rank].BestAv(cities)<< endl;
        }
    }
    if(rank == 0)
    {
        pop[rank].Chr[0].Print_cities(cities, path);
    }

    best[rank].close();
    av.close();
    path.close();
    MPI_Finalize();

    return 0;
}

void Initialisation(ofstream& best, int rank)
{
    best.open("path_best_mpi_core"+to_string(rank)+".dat0");
    av.open("path_average_mpi.dat0");
    path.open("best_final_path_mpi.dat0");

    // ====== INIZIALIZZAZIONE RND ======
    ifstream Primes, Seed, City;
    int p1, p2;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    // ====== LETTURA CITTA DA FILE ======
    City.open("American_capitals.dat");

    for(int i = 0; i < Ncities; i++)
    {   
        S_LEN += i+1;                             //per il check
        N_LEN += (i+1)*(i+1);
        double a, b;
        City >> a >> b;
        point c(a, b);
        cities.push_back(c);
    }
    City.close();
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
    for(int i = 0; i < Ncities; i++)
        M_Path[i]=i+1;
    for(int i = Ncities-1; i > 0; i--)      //algoritmo Fisher-Yates, l'altro non era più prestazionale
    {
        int j = (int)(rnd.Rannyu()*(i-1)+1);
        Swap(M_Path[i], M_Path[j]);
    }

    if(Check()==false){throw runtime_error("Fill Misfunction");}
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

void Population :: Birth(Random &rnd, int n)
{
    Npop = n;
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
        len += pow(cities[M_Path[i]-1].x()-cities[M_Path[i+1]-1].x(), 2) + pow(cities[M_Path[i]-1].y()-cities[M_Path[i+1]-1].y(),2);
        len += pow(cities[M_Path[Ncities-1]-1].x()-cities[M_Path[0]-1].x(), 2) + pow(cities[M_Path[Ncities-1]-1].y()-cities[M_Path[0]-1].y(), 2);
    }
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

void Population :: Mutate(Random &rnd, vector<point>& cities)
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
        Chr_appo[p].pair_permutation(); //pair permutation
        if(Chr_appo[p].Check()==false){throw runtime_error("Permutation Misfunction");}
        double bolz = exp((Chr[p].L2(cities)-Chr_appo[p].L2(cities))*Chr[p].beta);
        if(prob <= bolz) 
        {
            Chr[p]=Chr_appo[p];
            SAaccepted ++;
        }      

        prob = rnd.Rannyu();
        Chr_appo[p].n_shift();          //shift of n positions of m continguous cities  
        if(Chr_appo[p].Check()==false){throw runtime_error("Shift Misfunction");}
        bolz = exp((Chr[p].L2(cities)-Chr_appo[p].L2(cities))*Chr[p].beta);
        if(prob <= bolz)                
        {
            Chr[p]=Chr_appo[p];
            SAaccepted ++;
        }   

        prob = rnd.Rannyu();
        Chr_appo[p].m_permutation();    //permutation among m contiguous cities with other m different cities   
        if(Chr_appo[p].Check()==false){throw runtime_error("PermutaM Misfunction");}
        bolz = exp((Chr[p].L2(cities)-Chr_appo[p].L2(cities))*Chr[p].beta);
        if(prob <= bolz)                
        {
            Chr[p]=Chr_appo[p];
            SAaccepted ++;
        }   
        

        prob = rnd.Rannyu();
        Chr_appo[p].inversion();        //inversion of m cities
        if(Chr_appo[p].Check()==false){throw runtime_error("Inversion Misfunction");}
        bolz = exp((Chr[p].L2(cities)-Chr_appo[p].L2(cities))*Chr[p].beta);
        if(prob <= bolz)                
        {
            Chr[p]=Chr_appo[p];
            SAaccepted ++;
        }   
        SAattempted += 4;
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

Path& Path :: operator= (const Path& pth){  //scambia le stringhe di città, ma non cambia le temperature
    for(int i = 0; i < Ncities; i++)        //posso utilizzare la Template Swap per effettuare gli scambi nel programma
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

void SAReset (void) //Reset SA
{
  SAattempted = 0;
  SAaccepted = 0;
}

void Adjacent_Exchange(Path a, Path b)
{
    double bolz = exp(-(b.beta-a.beta)*(a.L2(cities)-b.L2(cities)));
    if(rnd.Rannyu() < bolz)
    {
        Swap(a, b);     //effettua lo scambio
        SAaccepted ++;  //aggiorna i contatori
    }
    SAattempted++;      //aggiorna i contatori
}

void Adjacent_Exchange_MPI(int &rank, int &size, Population pop, MPI_Status stat)
{
        int migrator[Ncities];  //vettore di appoggio, non possiamo spedire direttamente i Path
        for(int j = 0; j < Ncities; j++)
            migrator[j]=pop.Chr[nblk-1].M_Path[j];
        if(rank < size - 1)     //il rank spedisce il suo ultimo elemento al rank successivo
            MPI_Send(migrator, Ncities, MPI_INTEGER, rank + 1, 1, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        if(rank > 0)
        {                       //riceve l'elemento spedito dal rank che lo precede
            MPI_Recv(migrator, Ncities, MPI_INTEGER, rank - 1, 1, MPI_COMM_WORLD, &stat);
            Path Chr_neigh;     //Path di appoggio per il calcolo di L2
            for(int j = 0; j < Ncities; j++)
                Chr_neigh.M_Path[j]=migrator[j];
            double r = rnd.Rannyu();
            double bolz = exp(-(pop.Chr[0].L2(cities)+Chr_neigh.L2(cities))*(Chr_neigh.beta-pop.Chr[0].beta));
                                //valuta se effettuare lo scambio con probabilità in Lec08, p.30
            if(r <= bolz)       //in caso spedisce il suo primo elemento al rank precedente   
            {
            for(int j = 0; j < Ncities; j++)
                migrator[j]=pop.Chr[0].M_Path[j];
            MPI_Send(migrator, Ncities, MPI_INTEGER, rank - 1, 1, MPI_COMM_WORLD);
            pop.Chr[0] = Chr_neigh;
            SAaccepted ++;
            } else {            //altrimenti rispedisce al mittente
                MPI_Send(migrator, Ncities, MPI_INTEGER, rank - 1, 1, MPI_COMM_WORLD);
            }
            SAattempted ++;
            if(pop.Chr[0].Check()==false){throw runtime_error("Adjacent Exchange Misfunction: Send");}
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank < size - 1)     //il mittente riceve il risultato della scelta
        {
            MPI_Recv(migrator, Ncities, MPI_INTEGER, rank + 1, 1, MPI_COMM_WORLD, &stat);
            for(int j = 0; j < Ncities; j++)
                pop.Chr[nblk - 1].M_Path[j]=migrator[j];
                if(pop.Chr[nblk - 1].Check()==false){throw runtime_error("Adjacent Exchange Misfunction: Receive");}
        }
        MPI_Barrier(MPI_COMM_WORLD);
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