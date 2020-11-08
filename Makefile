CC     = gcc
CFLAGS = -O0 -g 

EASELDIR=/Users/dkaxiras/Desktop/lab/null_model/easel
all:   count_occurences

count_occurences: count_occurences.c Makefile
	${CC} ${CFLAGS} -o count_occurences -L ${EASELDIR}  -I ${EASELDIR}  count_occurences.c -leasel -lm -lpthread


clean:
	-rm *.o *~
	-rm count_occurrences
