FLAGS=-arch=compute_20 -code=sm_20
Warnings=-g -G

all: force.o euler.o
	nvcc main.cu force.o euler.o $(FLAGS) $(Warnings) -o test

force.o:
	nvcc force.cu $(FLAGS) $(Warnings) -c

euler.o:
	nvcc euler.cu $(FLAGS) $(Warnings) -c

clean:
	rm -rf *o test
