CC = g++
CFLAGS = -Wall -O3 --std=c++11

all : run1 run2

run2 : main2.exe
	./main2.exe
run1 : main.exe
	./main.exe

main2.exe : main2.o random.o VectorOperations.h
	$(CC) random.o main2.o -o main2.exe 
main2.o : main2.cpp 
	$(CC) -c main2.cpp -o main2.o $(CFLAGS) 
main.exe : main.o random.o IntegralMC.h
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main2.exe main.exe seed.out
