CC = gcc -g

.PHONY:	transform clean

trajectory:	main.o	coord.o
	$(CC)	-o transform main.o coord.o

main.o:	main.c
	$(CC)	-c main.c

coord.o:	coord.c
	$(CC)	-c coord.c

clean:
	-rm *.o
	-rm transform
