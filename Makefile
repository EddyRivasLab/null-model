CC     = gcc
CFLAGS = -O0 -g 

EASELDIR=/home/npcarter/hmmer/h4/lib/easel
all:   count_occurences

count_occurences: count_occurences.c Makefile
	${CC} ${CFLAGS} -o count_occurences -L ${EASELDIR}  -I ${EASELDIR}  count_occurences.c -leasel -lm -lpthread


clean:
	-rm *.o *~
	-rm count_occurrences
