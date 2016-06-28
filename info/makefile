EXECS=info
MPICC?=mpicc

all: ${EXECS}

info: info.c
	${MPICC} -o info info.c -std=c99

clean:
	rm -f ${EXECS}
