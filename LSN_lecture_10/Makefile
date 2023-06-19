CC = g++
CFLAGS = -Wall -O3 --std=c++11

run_mpi_6: TSP_MPI.exe
	mpiexec -np 6 TSP_MPI.exe

run: TSP.exe
	./TSP.exe

TSP_MPI.exe : TSP_MPI.cpp random.o VectorOperations.h TSP_MPI.h
	mpicxx TSP_MPI.cpp random.o -o TSP_MPI.exe $(CFLAGS)
TSP.exe : TSP.o random.o
	$(CC) random.o TSP.o -o TSP.exe
TSP.o : TSP.cpp TSP.h
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o TSP.exe TSP_MPI.exe

clean_dat0 : 
	rm *.dat0