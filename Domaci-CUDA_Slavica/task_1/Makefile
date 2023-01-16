# C compiler
CC = nvcc
CC_FLAGS = -O3
SEQ_CC = gcc

prime: 
	$(CC) $(CC_FLAGS) prime.cu -o prime

sequential:
	$(SEQ_CC) $(CC_FLAGS) prime.c -o prime_seq

clean:
	rm -f prime
	rm -f prime_seq
