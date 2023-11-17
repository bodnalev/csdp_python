#./Makefile
export CFLAGS=-m64 -march=native -mtune=native -Ofast -fopenmp -ansi -Wall -DBIT64 -DUSEOPENMP -DSETNUMTHREADS -DUSESIGTERM -DUSEGETTIME -I../include
export LIBS=-static -L../lib -lsdp -llapack -lblas -lm
all:
	cd lib; make libsdp.so
clean:
	cd lib; make clean