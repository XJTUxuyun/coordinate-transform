CC = gcc -g

.PHONY:	trajectory clean

trajectory:	main.o	coord.o
	$(CC)	-o trajectory main.o coord.o

main.o:	main.c
	$(CC)	-c main.c

coord.o:	coord.c
	$(CC)	-c coord.c

clean:
	-rm *.o
	-rm trajectory
