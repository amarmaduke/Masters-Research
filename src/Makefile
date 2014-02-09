FLAGS=-arch=compute_20 -code=sm_20
Warnings=-g -G
Links=-lcublas
C11=--compiler-options -std=gnu++0x

all: force.o euler.o dorpi.o utility.o
	nvcc main.cu force.o euler.o dorpi.o utility.o $(FLAGS) $(Warnings) $(Links) -o all.out

test: ./build/utility.o
	nvcc ./utility/test.cu $(FLAGS) $(Warnings) $(Links) -c -o ./build/utility_test.o
	nvcc ./test/test.cu ./build/utility.o ./build/utility_test.o $(FLAGS) $(Warnings) $(Links) -o test.out

force.o: utility.o
	nvcc force.cu $(FLAGS) $(Warnings) -c

euler.o: utility.o
	nvcc euler.cu $(FLAGS) $(Warnings) -c

dorpi.o: utility.o
	nvcc ./solver/dorpi.cu $(FLAGS) $(Warnings) -c -o ./build/dorpi.o

utility.o:
	nvcc ./utility/utility.cu $(FLAGS) $(Warnings) -c -o ./build/utility.o

clean:
	rm -rf ./build/*o
