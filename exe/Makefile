# makefile for petsc project
CC = /home/kongkong/Software/mpich-4.1.2/bin/mpicc
SRC = ../src
#CFLAGS = -c -std=gnu++20
CFLAGS = -c

CMD = -Wall -Wwrite-strings \
-Wno-strict-aliasing -Wno-unknown-pragmas \
-Wno-lto-type-mismatch -Wno-psabi -fstack-protector \
-fvisibility=hidden -g -O0

INCL = -I/home/kongkong/Software/petsc-3.19.2/include \
-I/home/kongkong/Software/mpich-4.1.2/include

LIB_PATH = -L/home/kongkong/Software/petsc-3.19.2/lib \
-Wl,-rpath,/usr/local/lib64 -L/usr/local/lib64 \
-L/home/kongkong/Software/mpich-4.1.2/lib \
-L/usr/lib/gcc/x86_64-redhat-linux/13 \
-Wl,-rpath,/home/kongkong/Software/mpich-4.1.2/lib \
-L/home/kongkong/LibraryObject/LIB_AGMG_PAR
#-L/home/kongkong/LibraryObject/LIB_AGMG

LIB_FLAG = -lAGMG_PAR -lpetsc -llapack -lblas -lm -lmpifort -lmpi \
-lgfortran -lgcc_s -lquadmath -lstdc++

LIB = $(LIB_PATH) $(LIB_FLAG)

OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: app_agmg

app_agmg: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_agmg $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm app_agmg $(OBJECT)
