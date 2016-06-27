EXECS=dot_product
MPICC?=mpicc

all: ${EXECS}

dot_product: dot_product.c
	${MPICC} -o dot_product dot_product.c

clean:
	rm -f ${EXECS}
