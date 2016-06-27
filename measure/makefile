EXECS=measure
MPICC?=mpicc

all: ${EXECS}

measure: measure.c
	${MPICC} -o measure measure.c

clean:
	rm -f ${EXECS}
