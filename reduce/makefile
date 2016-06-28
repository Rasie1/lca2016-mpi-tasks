EXECS=reduce
MPICC?=mpicc

all: ${EXECS}

reduce: reduce.c
	${MPICC} -o reduce reduce.c -std=c99

clean:
	rm -f ${EXECS}
