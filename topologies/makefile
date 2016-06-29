EXECS=topologies
MPICC?=mpicc

all: ${EXECS}

topologies: topologies.c
	${MPICC} -o topologies topologies.c -std=c99

clean:
	rm -f ${EXECS}
