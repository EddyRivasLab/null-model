CC     = gcc
CFLAGS = -O0 -g

EASELDIR=/Users/dkaxiras/Desktop/thesis/easel
all:   calculate_probabilities_with_output score_sequences

calculate_probabilities_with_output: calculate_probabilities_with_output.c Makefile
	${CC} ${CFLAGS} -o calculate_probabilities_with_output -L ${EASELDIR}  -I ${EASELDIR}  calculate_probabilities_with_output.c -leasel -lm -lpthread

score_sequences: score_sequences.c Makefile
	${CC} ${CFLAGS} -o score_sequences -L ${EASELDIR}  -I ${EASELDIR}  score_sequences.c -leasel -lm -lpthread

clean:
	-rm *.o *~
	-rm calculate_probabilities_with_output
	-rm score_sequences
