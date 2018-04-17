MATLAB_DIST_PATH="/media/aslab/data/matlab"
CC=$(shell which gcc)
BIN_FILES=$(shell find bin -type f)
BUILD_FILES=$(shell find build -type f)

all: buildfolders native mex

native: bin/testopencl LibSPH

mex: bin/mex/init.mexa64 bin/mex/callkernel.mexa64 bin/mex/computebins.mexa64 bin/mex/syncalldevicetohost.mexa64 bin/mex/syncallhosttodevice.mexa64 bin/mex/testinit.mexa64 bin/mex/querystate.mexa64 bin/mex/simstep.mexa64 bin/mex/setstate.mexa64 LibSPH_mex bin/mex/init_wrap.m bin/mex/testinit_wrap.m bin/mex/testinit_fluid.mexa64 bin/mex/simstep_fluid.mexa64 bin/mex/testinit_fluid_wrap.m

bin/test_conf_reorder: build/test_conf_reorder.o bin/libsph.so
	$(CC) -Wl,-rpath,'$$ORIGIN' -fPIC $< -o $@ -lsph -Lbin -lm

buildfolders:
	mkdir -p bin/mex
	mkdir -p build/3rdparty
	mkdir -p build/opencl
	mkdir -p build/mex/opencl

clean:
ifneq ($(strip $(BIN_FILES)),)
	rm $(BIN_FILES)
endif
ifneq ($(strip $(BUILD_FILES)),)
	rm $(BUILD_FILES)
endif

bin/testopencl: build/testopencl.o bin/libsph.so
	$(CC) -Wl,-rpath,'$$ORIGIN' -fPIC $< -o $@ -lsph -Lbin -lm -lOpenCL
	@echo "$$(tput bold)$$(tput setaf 2)Built $@$$(tput sgr0)"

LibSPH: bin/libsph.so

bin/libsph.so: build/particle_system.o build/opencl/particle_system_host.o build/opencl/platforminfo.o build/opencl/clerror.o build/note.o build/config.o build/build_psdata.o build/stringly.o build/3rdparty/whereami.o
	$(CC) -shared -fPIC $^ -o $@ -lOpenCL
	@echo "$$(tput bold)$$(tput setaf 2)Built $@$$(tput sgr0)"

build/%.o: src/%.c
	$(CC) -fPIC -c -I/media/aslab/data/matlab/extern/include $< -o $@

bin/mex/%.m: src/mex/%.m
	cp $< $@

LibSPH_mex: bin/mex/libsph_mex.so

bin/mex/libsph_mex.so: build/mex/particle_system_mex.o build/mex/opencl/particle_system_host_mex.o build/mex/opencl/platforminfo_mex.o build/mex/note_mex.o build/mex/config_mex.o build/opencl/clerror.o build/3rdparty/whereami.o build/mex/stringly_mex.o build/mex/build_psdata_mex.o
	$(CC) -shared -fPIC -I/media/aslab/data/matlab/extern/include $^ -o $@ -lOpenCL
	@echo "$$(tput bold)$$(tput setaf 2)Built $@$$(tput sgr0)"

bin/mex/%.mexa64: build/mex/%.o bin/mex/libsph_mex.so
	@$(CC) -pthread -Wl,--no-undefined -Wl,-rpath-link,$(MATLAB_DIST_PATH)/bin/glnxa64 \
		-shared  -O -Wl,--version-script,"$(MATLAB_DIST_PATH)/extern/lib/glnxa64/mexFunction.map" \
		$<  -Wl,-rpath,'$$ORIGIN' -lsph_mex   -Lbin/mex   -L"$(MATLAB_DIST_PATH)/bin/glnxa64" \
		-lmx -lmex -lmat -lm -lstdc++ -o $@
	@echo "$$(tput bold)$$(tput setaf 2)Built MEX file $@$$(tput sgr0)"

build/mex/%.o: src/mex/%.c
	$(MATLAB_DIST_PATH)/bin/mex -O CFLAGS="$$CFLAGS -fPIC -std=c99 -I/media/aslab/data/matlab/extern/include" -c $< -outdir build/mex

build/mex/%_mex.o: src/%.c
	$(MATLAB_DIST_PATH)/bin/mex -O CFLAGS="$$CFLAGS -fPIC -std=c99" -c $< -outdir build/mex/$(*D)
	mv build/mex/$*.o $@
