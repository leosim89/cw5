#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <vector>
#include <cstdio>
#include <iostream>
#include <string>
#include <cstdint>
#include <time.h>

#include <alloca.h>
#define __CL_ENABLE_EXCEPTIONS 
#include "CL/cl.hpp"

// includes for kernel text code
#include <fstream>
#include <streambuf>

#if !(defined(_WIN32) || defined(_WIN64))
#include <unistd.h>
void set_binary_io()
{}
#else
// http://stackoverflow.com/questions/341817/is-there-a-replacement-for-unistd-h-for-windows-visual-c
// http://stackoverflow.com/questions/13198627/using-file-descriptors-in-visual-studio-2010-and-windows
// Note: I could have just included <io.h> and msvc would whinge mightily, but carry on

#include <io.h>
#include <fcntl.h>

#define read _read
#define write _write
#define STDIN_FILENO 0
#define STDOUT_FILENO 1

void set_binary_io()
{
	_setmode(_fileno(stdin), _O_BINARY);
	_setmode(_fileno(stdout), _O_BINARY);
}
#endif


////////////////////////////////////////////
// Routines for bringing in binary images

/*! Reverse the orders of bits if necessary
\note This is laborious and a bit pointless. I'm sure it could be removed, or at least moved...
*/
uint64_t shuffle64(unsigned bits, uint64_t x)
{
	if (bits == 1){
		x = ((x & 0x0101010101010101ull) << 7)
			| ((x & 0x0202020202020202ull) << 5)
			| ((x & 0x0404040404040404ull) << 3)
			| ((x & 0x0808080808080808ull) << 1)
			| ((x & 0x1010101010101010ull) >> 1)
			| ((x & 0x2020202020202020ull) >> 3)
			| ((x & 0x4040404040404040ull) >> 5)
			| ((x & 0x8080808080808080ull) >> 7);
	}
	else if (bits == 2){
		x = ((x & 0x0303030303030303ull) << 6)
			| ((x & 0x0c0c0c0c0c0c0c0cull) << 2)
			| ((x & 0x3030303030303030ull) >> 2)
			| ((x & 0xc0c0c0c0c0c0c0c0ull) >> 6);
	}
	else if (bits == 4){
		x = ((x & 0x0f0f0f0f0f0f0f0full) << 4)
			| ((x & 0xf0f0f0f0f0f0f0f0ull) >> 4);
	}
	return x;
}

/*! Take data packed into incoming format, and exand to one integer per pixel */
void unpack_blob(unsigned w, unsigned h, unsigned bits, const uint64_t *pRaw, uint32_t *pUnpacked)
{
	uint64_t buffer = 0;
	unsigned bufferedBits = 0;

	const uint64_t MASK = 0xFFFFFFFFFFFFFFFFULL >> (64 - bits);

	for (unsigned i = 0; i<w*h; i++){
		if (bufferedBits == 0){
			buffer = shuffle64(bits, *pRaw++);
			bufferedBits = 64;
		}

		pUnpacked[i] = uint32_t(buffer&MASK);
		buffer = buffer >> bits;
		bufferedBits -= bits;
	}

	assert(bufferedBits == 0);
}

/*! Go back from one integer per pixel to packed format for output. */
void pack_blob(unsigned w, unsigned h, unsigned bits, const uint32_t *pUnpacked, uint64_t *pRaw)
{
	uint64_t buffer = 0;
	unsigned bufferedBits = 0;

	const uint64_t MASK = 0xFFFFFFFFFFFFFFFFULL >> (64 - bits);

	for (unsigned i = 0; i<w*h; i++){
		buffer = buffer | (uint64_t(pUnpacked[i] & MASK) << bufferedBits);
		bufferedBits += bits;

		if (bufferedBits == 64){
			*pRaw++ = shuffle64(bits, buffer);
			buffer = 0;
			bufferedBits = 0;
		}
	}

	assert(bufferedBits == 0);
}

bool read_blob(int fd, uint64_t cbBlob, void *pBlob)
{
	uint8_t *pBytes = (uint8_t*)pBlob;

	uint64_t done = 0;
	while (done<cbBlob){
		int todo = (int)std::min(uint64_t(1) << 30, cbBlob - done);

		int got = read(fd, pBytes + done, todo);
		if (got == 0 && done == 0)
			return false;	// end of file
		if (got <= 0){
			throw std::invalid_argument("Read failure.");
		}
		done += got;
	}

	return true;
}

void write_blob(int fd, uint64_t cbBlob, const void *pBlob)
{
	const uint8_t *pBytes = (const uint8_t*)pBlob;

	uint64_t done = 0;
	while (done<cbBlob){
		int todo = (int)std::min(uint64_t(1) << 30, cbBlob - done);

		int got = write(fd, pBytes + done, todo);
		if (got <= 0)
			throw std::invalid_argument("Write failure.");
		done += got;
	}
}

///////////////////////////////////////////////////////////////////
// Basic image processing primitives

uint32_t vmin(uint32_t a, uint32_t b)
{
	return std::min(a, b);
}

uint32_t vmin(uint32_t a, uint32_t b, uint32_t c)
{
	return std::min(a, std::min(b, c));
}

uint32_t vmin(uint32_t a, uint32_t b, uint32_t c, uint32_t d)
{
	return std::min(std::min(a, d), std::min(b, c));
}

uint32_t vmin(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e)
{
	return std::min(e, std::min(std::min(a, d), std::min(b, c)));
}

void erode_pixel(
	unsigned l,
	unsigned x,
	unsigned levels,
	unsigned h,
	const uint32_t *input,
	const uint32_t *right,
	const uint32_t *left,
	uint32_t *output) {

	int index = l*(h + 1) + x;

	if (input[(l + 1) * h + l]){

		if (!left[(l + 1) * h + l]) {
			if (x == 0) 
				output[index + (h+1)] = vmin(input[index], input[index + 1], right[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmin(input[index], input[index - 1], right[index]);
				output[index + (h + 2)] = 1;
			}
			else 
				output[index + (h + 1)] = vmin(input[index], input[index + 1], input[index - 1], right[index]);
		}

		else if (!right[(l + 1) * h + l]) {
			if (x == 0)
				output[index + (h + 1)] = vmin(input[index], input[index + 1], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmin(input[index], input[index - 1], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmin(input[index], input[index + 1], input[index - 1], left[index]);
		}

		else {
			if (x == 0)
				output[index + (h + 1)] = vmin(input[index], input[index + 1], right[index], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmin(input[index], input[index - 1], right[index], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmin(input[index], input[index + 1], input[index - 1], right[index], left[index]);
		}

	} 
	else {
		output[index + (h + 2)] = 0;
	}


}


uint32_t vmax(uint32_t a, uint32_t b)
{
	return std::max(a, b);
}

uint32_t vmax(uint32_t a, uint32_t b, uint32_t c)
{
	return std::max(a, std::max(b, c));
}

uint32_t vmax(uint32_t a, uint32_t b, uint32_t c, uint32_t d)
{
	return std::max(std::max(a, d), std::max(b, c));
}

uint32_t vmax(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e)
{
	return std::max(e, std::max(std::max(a, d), std::max(b, c)));
}

void dilate_pixel(
	unsigned l,
	unsigned x,
	unsigned levels,
	unsigned h,
	const uint32_t *input,
	const uint32_t *right,
	const uint32_t *left,
	uint32_t *output) {

	int index = l*(h + 1) + x;

	if (input[(l + 1) * h + l]){

		if (!left[(l + 1) * h + l]) {
			if (x == 0)
				output[index + (h + 1)] = vmax(input[index], input[index + 1], right[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmax(input[index], input[index - 1], right[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmax(input[index], input[index + 1], input[index - 1], right[index]);
		}

		else if (!right[(l + 1) * h + l]) {
			if (x == 0)
				output[index + (h + 1)] = vmax(input[index], input[index + 1], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmax(input[index], input[index - 1], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmax(input[index], input[index + 1], input[index - 1], left[index]);
		}

		else {
			if (x == 0)
				output[index + (h + 1)] = vmax(input[index], input[index + 1], right[index], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmax(input[index], input[index - 1], right[index], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmax(input[index], input[index + 1], input[index - 1], right[index], left[index]);
		}

	}
	else {
		output[index + (h + 2)] = 0;
	}

}

///////////////////////////////////////////////////////////////////
// Composite image processing

// You may want to play with this to check you understand what is going on
void invert(unsigned w, unsigned h, unsigned bits, std::vector<uint32_t> &pixels)
{
	uint32_t mask = 0xFFFFFFFFul >> bits;

	for (unsigned i = 0; i<w*h; i++){
		pixels[i] = mask - pixels[i];
	}
}

void kernel_lx (
	unsigned l,
	unsigned x,
	int levels,
	unsigned h,
	const uint32_t *input,
	const uint32_t *right,
	const uint32_t *left,
	uint32_t *output)
{
	unsigned absl = std::abs(levels);
	if (levels<0) {
		erode_pixel(l, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
		dilate_pixel(l + absl, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
	} else {
		dilate_pixel(l, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
		erode_pixel(l + absl, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
	}
}

//function to load kernel code
std::string LoadSource(const char *fileName)
{
    std::string baseDir="src";
    if(getenv("HPCE_CL_SRC_DIR")){
        baseDir=getenv("HPCE_CL_SRC_DIR");
    }
	
	std::string fullName=baseDir+"/"+fileName;
	 
	std::ifstream src(fullName, std::ios::in | std::ios::binary);
	if(!src.is_open())
		throw std::runtime_error("LoadSource : Couldn't load cl file from '"+fullName+"'.");
		
	return std::string(
        (std::istreambuf_iterator<char>(src)), // Node the extra brackets.
        std::istreambuf_iterator<char>()
    );
}

int main(int argc, char *argv[])
{
	try{
		if (argc<3){
			fprintf(stderr, "Usage: process width height [bits] [levels]\n");
			fprintf(stderr, "   bits=8 by default\n");
			fprintf(stderr, "   levels=1 by default\n");
			exit(1);
		}

		unsigned w = atoi(argv[1]);
		unsigned h = atoi(argv[2]);

		unsigned bits = 8;
		if (argc>3){
			bits = atoi(argv[3]);
		}

		if (bits>32)
			throw std::invalid_argument("Bits must be <= 32.");

		unsigned tmp = bits;
		while (tmp != 1){
			tmp >>= 1;
			if (tmp == 0)
				throw std::invalid_argument("Bits must be a binary power.");
		}

		if (((w*bits) % 64) != 0){
			throw std::invalid_argument(" width*bits must be divisible by 64.");
		}

		int levels = 1;
		if (argc>4){
			levels = atoi(argv[4]);
		}
		unsigned absl = std::abs(levels);

		fprintf(stderr, "Processing %d x %d image with %d bits per pixel.\n", w, h, bits);

		/************************OpenCL******************************/
			//Choosing a platform
		std::vector<cl::Platform> platforms;

		cl::Platform::get(&platforms);
		if(platforms.size()==0)
			throw std::runtime_error("No OpenCL platforms found.");
		
		std::cerr<<"Found "<<platforms.size()<<" platforms\n";
		for(unsigned i=0;i<platforms.size();i++){
			std::string vendor=platforms[0].getInfo<CL_PLATFORM_VENDOR>();
			std::cerr<<"  Platform "<<i<<" : "<<vendor<<"\n";
		}
	
		int selectedPlatform=0;
		if(getenv("HPCE_SELECT_PLATFORM")){
			selectedPlatform=atoi(getenv("HPCE_SELECT_PLATFORM"));
		}
		std::cerr<<"Choosing platform "<<selectedPlatform<<"\n";
		cl::Platform platform=platforms.at(selectedPlatform);  
	
		//Choosing a device
		std::vector<cl::Device> devices;
		platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);  
		if(devices.size()==0){
			throw std::runtime_error("No opencl devices found.\n");
		}

		std::cerr<<"Found "<<devices.size()<<" devices\n";
		for(unsigned i=0;i<devices.size();i++){
			std::string name=devices[i].getInfo<CL_DEVICE_NAME>();
			std::cerr<<"  Device "<<i<<" : "<<name<<"\n";
		}
	
		int selectedDevice=0;
		if(getenv("HPCE_SELECT_DEVICE")){
			selectedDevice=atoi(getenv("HPCE_SELECT_DEVICE"));
		}
		std::cerr<<"Choosing device "<<selectedDevice<<"\n";
		cl::Device device=devices.at(selectedDevice);
	
		cl::Context context(devices); 

		std::string kernelSource=LoadSource("process_v4_kernel.cl");

		//compiling sources into context
		cl::Program::Sources sources;   // A vector of (data,length) pairs
		sources.push_back(std::make_pair(kernelSource.c_str(), kernelSource.size()+1)); // push on our single string
	
		cl::Program program(context, sources);
		//catching errors thrown to view errors in cmd
		try{
			program.build(devices);
		}catch(...){
			for(unsigned i=0;i<devices.size();i++){
			std::cerr<<"Log for device "<<devices[i].getInfo<CL_DEVICE_NAME>()<<":\n\n";
			std::cerr<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[i])<<"\n\n";
			}
			throw;
		}

		//allocating buffers
		size_t cbBuffer=4*(((absl * 2) + 1) * (h + 1));
		cl::Buffer buffOutput(context, CL_MEM_WRITE_ONLY, cbBuffer);
		cl::Buffer buffCenter(context, CL_MEM_READ_ONLY, cbBuffer);
		cl::Buffer buffRight(context, CL_MEM_READ_ONLY, cbBuffer);
		cl::Buffer buffLeft(context, CL_MEM_READ_ONLY, cbBuffer);

		/*******************************************************************/
		uint64_t cbRaw = uint64_t(1)*h*bits / 8;
		std::vector<uint64_t> raw(size_t(cbRaw / 8));
		std::vector<uint32_t> output(((absl * 2) + 1) * (h + 1));
		std::vector<uint32_t> center(((absl * 2) + 1) * (h + 1));
		std::vector<uint32_t> right(((absl * 2) + 1) * (h + 1));
		std::vector<uint32_t> left(((absl * 2) + 1) * (h + 1));
		for (unsigned i = 0; i < (((absl * 2) + 1) * (h + 1)); i++){
			output[i] = 0;
			center[i] = 0;
			right[i] = 0;
			left[i] = 0;
		}

		set_binary_io();

		clock_t *t;
		t = new clock_t[w];

		/************KERNEL**********************/
		//finding the compiled kernel and creating an instance of it
		cl::Kernel kernel(program, "kernel_lx");

		//Setting the Kernel Params  --> can bring outside for loop?
		kernel.setArg(0, levels);
		kernel.setArg(1, h);
		kernel.setArg(2, buffCenter);
		kernel.setArg(3, buffRight);
		kernel.setArg(4, buffLeft);
		kernel.setArg(5, buffOutput);
		//creating command queue for single device
		cl::CommandQueue queue(context, device);

		//Setting up the iteration space
		cl::NDRange offset(0, 0);               // Always start iterations at x=0, y=0
		cl::NDRange globalSize(absl, h);   // Global size must match the original loops
		cl::NDRange localSize=cl::NullRange;

		while (1) {

			t[0] = clock();

			right[h] = 1;
			if (!read_blob(STDIN_FILENO, cbRaw, &raw[0]))
				break;
			unpack_blob(1, h, bits, &raw[0], &right[0]);

			if (absl) {
				for (unsigned i = 1; i < (w + 4*absl-1); i++){ 

					// Managing the input and output buffers in the pipeline
					std::swap(left, center);
					std::swap(center, right);
					std::swap(output, right);

					if (i < w){
						t[i] = clock();
						right[h] = 1;
						if (!read_blob(STDIN_FILENO, cbRaw, &raw[0]))
							break;
						unpack_blob(1, h, bits, &raw[0], &right[0]);
					} else {
						right[h] = 0;
					}

					/**************KERNEL****************/
					queue.enqueueWriteBuffer(buffCenter, CL_TRUE, 0, cbBuffer, &center[0]);
					queue.enqueueWriteBuffer(buffRight, CL_TRUE, 0, cbBuffer, &right[0]);
					queue.enqueueWriteBuffer(buffLeft, CL_TRUE, 0, cbBuffer, &left[0]);

					queue.enqueueNDRangeKernel(kernel, offset, globalSize, localSize);
					queue.enqueueBarrier();
					queue.enqueueReadBuffer(buffOutput, CL_TRUE, 0, cbBuffer, &output[0]);

					// Capture output
					if (i >= (4*absl-1)){
						pack_blob(1, h, bits, &output[(absl*2) * (h+1)], &raw[0]);
						write_blob(STDOUT_FILENO, cbRaw, &raw[0]);
						t[i - 4 * absl + 1] = clock() - t[i - 4 * absl + 1];
					}

				}
			}
			else {
				pack_blob(1, h, bits, &right[0], &raw[0]);
				write_blob(STDOUT_FILENO, cbRaw, &raw[0]);
				for (unsigned i = 1; i < w; i++){ 
					if (!read_blob(STDIN_FILENO, cbRaw, &raw[0]))
						break;
					unpack_blob(1, h, bits, &raw[0], &right[0]);
					pack_blob(1, h, bits, &right[0], &raw[0]);
					write_blob(STDOUT_FILENO, cbRaw, &raw[0]);
				}
			}


			fprintf(stderr, "Cycle completed! \n");

			clock_t maxcl = t[0];
			for (unsigned int i = 1; i < w; i++) {
				if (t[i] > maxcl) maxcl = t[i];
			}

			fprintf(stderr, "Pixel Latency: %f \n", ((float)maxcl) / CLOCKS_PER_SEC);
		}
		return 0;
	}
	catch (std::exception &e){
		std::cerr << "Caught exception : " << e.what() << "\n";
		return 1;
	}
}