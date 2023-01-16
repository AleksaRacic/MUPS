# Makefile for the molecular dynamics code

#
# C compiler and options for GNU 
#
CC=/usr/bin/nvcc
LIB=-lm

#
# Object files
#
OBJ=    main.o \
	dfill.o \
	domove.o \
        dscal.o \
	fcc.o \
	forces.o \
	mkekin.o \
	mxwell.o \
	prnout.o \
	velavg.o

#
# Compile
#
md:
	cuda main.cu dfill.cu domove.cu dscal.cu fcc.cu forces.cu mkekin.cu mxwell.cu prnout.cu velavg.cu -o md -arch=sm_61 -lm

# md:	$(OBJ)
# 	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o md
