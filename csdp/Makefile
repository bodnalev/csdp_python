#./csdp/Makefile
export CFLAGS=-m64 -march=native -mtune=native -Ofast -fopenmp -ansi -Wall -DBIT64 -DUSEOPENMP -DSETNUMTHREADS -DUSESIGTERM -DUSEGETTIME -I../include -g
export LIBS=-static -L../src -lsdp -llapack -lblas -lm

all:
	cd src; make libsdp.so
	cd src; make libsdp.a

clean:
	cd src; make clean