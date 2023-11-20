#./Makefile
export CFLAGS=-m64 -march=native -mtune=native -Ofast -fopenmp -ansi -Wall -DBIT64 -DUSEOPENMP -DSETNUMTHREADS -DUSESIGTERM -DUSEGETTIME -I../include
export LIBS=-static -L../lib -lsdp -llapack -lblas -lm

all:
	cd lib; make libsdp.so
	cd lib; make libsdp.a
	cd test; make small_test
	cd test; make csdp

clean:
	cd lib; make clean
	cd test; make clean