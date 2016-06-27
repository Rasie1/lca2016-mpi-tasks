EXECS=big
MPICC?=mpicc

all: ${EXECS}

big: big.c
	${MPICC} -o big big.c

clean:
	rm -f ${EXECS}
