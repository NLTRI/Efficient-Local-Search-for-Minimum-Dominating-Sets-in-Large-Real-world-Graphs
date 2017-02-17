all: SAMDS

SAMDS: SAMDS.cpp
	g++ -std=gnu++0x -O3 -static SAMDS.cpp -o SAMDS

clean: rm -f *~
