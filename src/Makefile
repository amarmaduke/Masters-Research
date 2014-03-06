FLAGS=-arch=compute_20 -code=sm_20 -I/home/marmaduke/include/ -L/home/marmaduke/lib
BoostFix=-DBOOST_NOINLINE='__attribute__ ((noinline))'
Warnings=-g -G
Links=-lcublas
C11=--compiler-options -std=gnu++0x

odeint: force.o utility.o json.o json_checker.o
	nvcc ./simulation/main.cu ./build/force.o ./build/utility.o ./build/json.o ./build/json_checker.o $(BoostFix) $(FLAGS) $(Warnings) $(Links) -o odeint.out

home: force.o utility.o json.o json_checker.o
	nvcc ./simulation/main.cu ./build/force.o ./build/utility.o ./build/json.o ./build/json_checker.o -I ~/Research/Masters-Research/src/ $(FLAGS) $(Warnings) $(Links) -o odeint.out

all: force.o dorpi.o utility.o
	nvcc ./simulation/main.cu ./build/force.o ./build/dorpi.o ./build/utility.o $(FLAGS) $(Warnings) $(Links) -o all.out

test: utility.o
	nvcc ./utility/test.cu $(FLAGS) $(Warnings) $(Links) -c -o ./build/utility_test.o
	nvcc ./test/test.cu ./build/utility.o ./build/utility_test.o $(FLAGS) $(Warnings) $(Links) -o test.out

force.o: utility.o
	nvcc ./simulation/force.cu $(FLAGS) $(Warnings) -c -o ./build/force.o

dorpi.o: utility.o
	nvcc ./solver/dorpi.cu $(FLAGS) $(Warnings) -c -o ./build/dorpi.o

utility.o:
	nvcc ./utility/utility.cu $(FLAGS) $(Warnings) -c -o ./build/utility.o

json_checker.o:
	nvcc ./json/JSON_checker.cu $(FLAGS) $(Warnings) -c -o ./build/json_checker.o

json.o: json_checker.o
	nvcc ./json/json.cu $(FLAGS) $(Warnings) -c -o ./build/json.o

clean:
	rm -rf ./build/*o
	rm -rf ./*.out
