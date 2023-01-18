# C compiler
CC = nvcc
CC_FLAGS = -O3 -arch=sm_61
SEQ_CC = gcc

cuda: 
	$(CC) $(CC_FLAGS) feyman.cu -o feyman

sequential:
	$(SEQ_CC) -O3 feyman.c -o feyman_seq -lm

clean:
	rm -f feyman
	rm -f feyman_seq
