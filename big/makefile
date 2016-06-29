EXECS=big
MPICC?=mpicc

all: ${EXECS}

big: big.c
	${MPICC} -o big big.c -std=c99

clean:
	rm -f ${EXECS}
