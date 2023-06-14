CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: SA.exe
	./SA.exe
	
SA.exe : SA.o random.o
	$(CC) random.o SA.o -o SA.exe
SA.o : main.cpp
	$(CC) -c main.cpp -o SA.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SA.exe
clean_dat :
	rm *.dat
