CC     = gcc
CFLAGS = -O0 -g 

EASELDIR=/Users/dkaxiras/Desktop/thesis/easel
all:   count_occurences count_occurences_without_probabs

count_occurences: count_occurences.c Makefile
	${CC} ${CFLAGS} -o count_occurences -L ${EASELDIR}  -I ${EASELDIR}  count_occurences.c -leasel -lm -lpthread

count_occurences_without_probabs: count_occurences_without_probabs.c Makefile
	${CC} ${CFLAGS} -o count_occurences_without_probabs -L ${EASELDIR}  -I ${EASELDIR}  count_occurences_without_probabs.c -leasel -lm -lpthread

clean:
	-rm *.o *~
	-rm count_occurrences
	-rm count_occurences_without_probabs
