FLAGS=-arch=compute_20 -code=sm_20
Warnings=-g -G
Links=-lcublas

all: force.o euler.o dorpi.o
	nvcc main.cu force.o euler.o dorpi.o $(FLAGS) $(Warnings) $(Links) -o test

force.o:
	nvcc force.cu $(FLAGS) $(Warnings) -c

euler.o:
	nvcc euler.cu $(FLAGS) $(Warnings) -c

dorpi.o:
	nvcc dorpi.cu $(FLAGS) $(Warnings) -c

clean:
	rm -rf *o test
