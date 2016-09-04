MY_NAME=frqo
MY_CC=mpicc
MY_SYMBS=-DKEEP_FILES
MY_CFLAGS=-Wall -O3 -std=gnu99
MY_INCS=-I./ -I..
MY_LIBS=-L./
MY_LINKS=-lqform -loptarith -lgmp

all: myfunctions.c generator.c main.c 
	$(MY_CC) $(MY_CFLAGS) $(MY_SYMBS) $(MY_INCS) -o $(MY_NAME) $^ $(MY_LIBS) $(MY_LINKS)
lib:
	$(MY_CC) $(MY_CFLAGS) $(MY_SYMBS) $(MY_INCS) -c myfunctions.c -o myfunctions.o
	ar rcs lib$(MY_NAME).a myfunctions.o

clean:
	rm *.o
