Flags=-arch=compute_20 -code=sm_20 -I/home/marmaduke/include/ -L/home/marmaduke/lib
BoostFix=-DBOOST_NOINLINE='__attribute__ ((noinline))'
Debug=-g -G
Opt=-02 --compiler-options -fno-strict-aliasing
BuildDir=/home/marmaduke/Masters-Research/src/build/bepp
Dir=/home/marmaduke/Masters-Research/src/bepp/

debug: $(BuildDir)/force_debug.o $(BuildDir)/json.o $(BuildDir)/json_checker.o
	nvcc $(Dir)/routine.cu $^ $(BoostFix) $(Flags) $(Debug) -o ../bepp_debug.out

release: $(BuildDir)/json.o $(BuildDir)/json_checker.o
	nvcc $(Dir)/routine.cu $^ $(BoostFix) $(Flags) $(Opt) -o ../bepp.out

force: $(BuildDir)/force.o
$(BuildDir)/force.o:
	nvcc $(Dir)/force.cu $(BoostFix) $(Flags) $(Opt) -c -o $@

force_debug: $(BuildDir)/force_debug.o
$(BuildDir)/force_debug.o:
	nvcc $(Dir)/force.cu $(BoostFix) $(Flags) $(Warnings) -c -o $@

json: $(BuildDir)/json.o $(BuildDir)/json_checker.o
$(BuildDir)/json.o: $(BuildDir)/json_checker.o
	nvcc ../json/json.cu $(Flags) $(Debug) -c -o $@

$(BuildDir)/json_checker.o:
	nvcc ../json/JSON_checker.cu $(Flags) $(Debug) -c -o $@

setup:
	mkdir $(BuildDir)

clean:
	rm -rf $(BuildDir)/*o
	rm -rf ../bepp*.out
