CC     = gcc
CFLAGS = -O0 -g

EASELDIR=/Users/dkaxiras/Desktop/thesis/easel
all:   count_occurences count_occurences_with_output calculate_probabilities calculate_probabilities_with_output score_sequences

count_occurences: count_occurences.c Makefile
	${CC} ${CFLAGS} -o count_occurences -L ${EASELDIR}  -I ${EASELDIR}  count_occurences.c -leasel -lm -lpthread

count_occurences_with_output: count_occurences_with_output.c Makefile
	${CC} ${CFLAGS} -o count_occurences_with_output -L ${EASELDIR}  -I ${EASELDIR}  count_occurences_with_output.c -leasel -lm -lpthread

calculate_probabilities: calculate_probabilities.c Makefile
	${CC} ${CFLAGS} -o calculate_probabilities -L ${EASELDIR}  -I ${EASELDIR}  calculate_probabilities.c -leasel -lm -lpthread

calculate_probabilities_with_output: calculate_probabilities_with_output.c Makefile
	${CC} ${CFLAGS} -o calculate_probabilities_with_output -L ${EASELDIR}  -I ${EASELDIR}  calculate_probabilities_with_output.c -leasel -lm -lpthread

score_sequences: score_sequences.c Makefile
	${CC} ${CFLAGS} -o score_sequences -L ${EASELDIR}  -I ${EASELDIR}  score_sequences.c -leasel -lm -lpthread

clean:
	-rm *.o *~
	-rm count_occurences
	-rm count_occurences_with_output
	-rm calculate_probabilities
	-rm calculate_probabilities_with_output
	-rm score_sequences
