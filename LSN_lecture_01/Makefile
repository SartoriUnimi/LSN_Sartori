CC = g++
CFLAGS = -Wall -O3 --std=c++11 

all: main1.exe main2.exe main3.exe main4.exe main5.exe run1 run2 run3 run4 run5

run1:
	./main1.exe
run2:
	./main2.exe
run3:
	./main3.exe
run4:
	./main4.exe
run5:
	./main5.exe

main1.exe : main1.o random.o Functions.h VectorOperations.h
	$(CC) random.o main1.o -o main1.exe
main1.o : main1.cpp
	$(CC) -c main1.cpp -o main1.o $(CFLAGS)
main2.exe : main2.o random.o Functions.h VectorOperations.h
	$(CC) random.o main2.o -o main2.exe
main2.o : main2.cpp
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)
main3.exe : main3.o random.o Functions.h VectorOperations.h
	$(CC) random.o main3.o -o main3.exe
main3.o : main3.cpp
	$(CC) -c main3.cpp -o main3.o $(CFLAGS)
main4.exe : main4.o random.o Functions.h VectorOperations.h
	$(CC) random.o main4.o -o main4.exe
main4.o : main4.cpp
	$(CC) -c main4.cpp -o main4.o $(CFLAGS)
main5.exe : main5.o random.o Functions.h VectorOperations.h
	$(CC) random.o main5.o -o main5.exe
main5.o : main5.cpp
	$(CC) -c main5.cpp -o main5.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main1.exe main2.exe main3.exe main4.exe main5.exe seed.out
