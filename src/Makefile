FLAGS=-arch=compute_20 -code=sm_20 -I/home/marmaduke/include/ -L/home/marmaduke/lib
BoostFix=-DBOOST_NOINLINE='__attribute__ ((noinline))'
Warnings=-g -G
Links=-lcublas -lboost_thread -lboost_system
C11=--compiler-options -std=gnu++0x
BuildDir=/home/marmaduke/Masters-Research/src/build

odeint: $(BuildDir)/force.o $(BuildDir)/json.o $(BuildDir)/json_checker.o $(BuildDir)/main.o
	nvcc ./simulation/test.cpp $^ $(BoostFix) $(FLAGS) $(Warnings) $(Links) -o odeint.out

test2: jacobi json.o json_checker.o
	nvcc ./simulation/test.cpp ./build/jacobi.o ./build/json.o ./build/json_checker.o $(BoostFix) $(Flags) $(Warnings) -o test2.out

jacobi: utility.o
	nvcc ./simulation/jacobi.cu $(BoostFix) $(FLAGS) $(Warnings) -c -o ./build/jacobi.o

$(BuildDir)/main.o: $(BuildDir)/utility.o $(BuildDir)/json.o $(BuildDir)/force.o
	nvcc ./simulation/main.cu $(BoostFix) $(FLAGS) $(Warnings) -c -o $@

home: force.o utility.o json.o json_checker.o
	nvcc ./simulation/main.cu ./build/force.o ./build/utility.o ./build/json.o ./build/json_checker.o -I ~/Research/Masters-Research/src/ $(FLAGS) $(Warnings) $(Links) -o odeint.out

force: $(BuildDir)/force.o
$(BuildDir)/force.o: $(BuildDir)/utility.o
	nvcc ./simulation/force.cu $(FLAGS) $(Warnings) -c -o $@

$(BuildDir)/dorpi.o: $(BuildDir)/utility.o
	nvcc ./solver/dorpi.cu $(FLAGS) $(Warnings) -c -o $@

utility: $(BuildDir)/utility.o
$(BuildDir)/utility.o:
	nvcc ./utility/utility.cu $(FLAGS) $(Warnings) -c -o $@

json: $(BuildDir)/json.o $(BuildDir)/json_checker.o
$(BuildDir)/json.o: $(BuildDir)/json_checker.o
	nvcc ./json/json.cu $(FLAGS) $(Warnings) -c -o $@

$(BuildDir)/json_checker.o:
	nvcc ./json/JSON_checker.cu $(FLAGS) $(Warnings) -c -o $@

clean:
	rm -rf ./build/*o
	rm -rf ./*.out
