process_v0:
	g++ -std=c++11 -O2 src/process.cpp -o ./bin/process
	
process_v1_pipeline:
	g++ -std=c++11 -O2 src/process_v1_pipeline.cpp -o ./bin/process_v1_pipeline

process_v2_1D:
	g++ -std=c++11 -O2 src/process_v2_1D.cpp -o ./bin/process_v2_1D
	
process_v4_kernel:
	g++ -std=c++11 -O2 src/process_v4_kernel.cpp -I ./include/ -I ./opencl_sdk/include -L ./opencl_sdk/lib/cygwin/x86 -l OpenCL -o ./bin/process_v4_kernel

all: process_v0 process_v1_pipeline process_v2_1D process_v4_kernel