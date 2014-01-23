FLAGS=-arch=compute_20 -code=sm_20

all: force.o euler.o
	nvcc main.cu force.o euler.o $(FLAGS) -o test

force.o:
	nvcc force.cu $(FLAGS) -c

euler.o:
	nvcc euler.cu $(FLAGS) -c

clean:
	rm -rf *o test
