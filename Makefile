all: program

program: main.o particle.o
	g++ main.o particle.o -o program

main.o: main.cpp
	g++ -static -c main.cpp

particle.o: particle.cpp
	g++ -static -c particle.cpp

clean:
	rm -rf *.o program