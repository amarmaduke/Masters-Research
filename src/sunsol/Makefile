CUDA_FLAGS = -arch=compute_20 -code=sm_20 -I/home/marmaduke/include/ -L/home/marmaduke/lib
BoostFix=-DBOOST_NOINLINE='__attribute__ ((noinline))'
CUDA_Debug = -g -G
Debug = -g
CUDA_Opt = -O3 --compiler-options -fno-strict-aliasing
Opt = -03 -fno-strict-aliasing
BigBox_BuildDir=/home/marmaduke/Masters-Research/src/build/sunsol/
BigBox_Dir=/home/marmaduke/Masters-Research/src/bepp/
Home_BuildDir=/home/owner/Research/Masters-Research/src/build/sunsol/
Home_Dir=/home/owner/Research/Masters-Research/src/sunsol/
FLAGS = -I$(Home_Dir) -lsundials_cvode -lsundials_nvecserial -lm

debug_bigbox: $(BigBox_BuildDir)/force_debug.o $(BigBox_BuildDir)/json.o $(BigBox_BuildDir)/json_checker.o
	nvcc $(BigBox_Dir)/main.cpp $^ $(CUDA_FLAGS) $(CUDA_Debug) -o $(BigBox_Dir)/../sunsol_debug.out

release_bigbox: $(BigBox_BuildDir)/force.o $(BigBox_BuildDir)/json.o $(BigBox_BuildDir)/json_checker.o
	nvcc $(BigBox_Dir)/main.cpp $^ $(CUDA_FLAGS) $(CUDA_Opt) -o $(BigBox_Dir)/../sunsol.out

debug_home: $(Home_BuildDir)/force_debug.o $(Home_BuildDir)/json.o $(Home_BuildDir)/json_checker.o
	g++ $(Home_Dir)/main.cpp $^ $(Debug) $(FLAGS) -o $(Home_Dir)/../sunsol_debug.out

release_home: $(Home_BuildDir)/force.o $(Home_BuildDir)/json.o $(Home_BuildDir)/json_checker.o
	g++ $(Home_Dir)/main.cpp $^ $(Opt) $(FLAGS) -o $(Home_Dir)/../sunsol.out

force_bbox: $(BigBox_BuildDir)/force_debug.o
$(BigBox_BuildDir)/force.o:
	nvcc $(BigBox_Dir)/cuda/force.cu $(CUDA_FLAGS) $(CUDA_Opt) -c -o $@

force_bbox_debug: $(BigBox_BuildDir)/force_debug.o
$(BigBox_BuildDir)/force_debug.o:
	nvcc $(BigBox_Dir)/cuda/force.cu $(CUDA_FLAGS) $(CUDA_Debug) -c -o $@

force_home: $(Home_BuildDir)/force_debug.o
$(Home_BuildDir)/force.o:
	g++ $(Home_Dir)/serial/force.cpp $(Opt) $(FLAGS) -c -o $@

force_home_debug: $(Home_BuildDir)/force_debug.o
$(Home_BuildDir)/force_debug.o:
	g++ $(Home_Dir)/serial/force.cpp $(Debug) $(FLAGS) -c -o $@

json_bbox: $(BigBox_BuildDir)/json.o $(BigBox_BuildDir)/json_checker.o
$(BigBox_BuildDir)/json.o:
	g++ $(BigBox_Dir)/json/json.cpp -c -o $@

$(BigBox_BuildDir)/json_checker.o:
	g++ $(BigBox_Dir)/json/JSON_checker.cpp -c -o $@

json_home: $(Home_BuildDir)/json.o $(Home_BuildDir)/json_checker.o
$(Home_BuildDir)/json.o:
	g++ $(Home_Dir)/json/json.cpp -c -o $@

$(Home_BuildDir)/json_checker.o:
	g++ $(Home_Dir)/json/JSON_checker.cpp -c -o $@

clean_home:
	rm -rf $(Home_BuildDir)/*o
	rm -rf ../sunsol*.out

clean_bbox:
	rm -rf $(BigBox_BuildDir)/*o
	rm -rf ../sunsol*.out
