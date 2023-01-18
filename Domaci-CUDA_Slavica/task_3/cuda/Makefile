# Makefile for the molecular dynamics code

#
# Compile
#
md:	$(OBJ)
	nvcc main.cu dfill.cu domove.cu dscal.cu fcc.cu forces.cu mkekin.cu mxwell.cu prnout.cu velavg.cu -o md -arch=sm_61 -lm

clean:
	rm *.o md