#include "passPICuda.h"

#include "iNau.h"
#include "nau.h"
#include "nau/render/passFactory.h"
#include "nau/material/iTexture.h"


#include <glbinding/gl/gl.h>
#include <glbinding/Binding.h>

#ifdef WIN32
#include <Windows.h>
#endif

#include <cuda_gl_interop.h>
#include <cuda.h>


//using namespace gl;

//static nau::INau *NAU_INTERFACE;

static char className[] = "cudaPI";

#ifdef WIN32
#define EXPORT __declspec(dllexport)
#elif __linux__
#define EXPORT extern "C"
#endif

EXPORT 
void *
createPass(const char *s) {

	std::shared_ptr<PAssCudaPI> *p = new std::shared_ptr<PAssCudaPI>(new PAssCudaPI(s));
	return p;
}


EXPORT
void 
init(void *nauInst) {

	INau::SetInterface((nau::INau *)nauInst);
	nau::Nau::SetInstance((nau::Nau *)nauInst);
	glbinding::Binding::initialize(false);
}


EXPORT
char *
getClassName() {

	return className;
}


using namespace nau::geometry;
using namespace nau::math;
using namespace nau::render;
using namespace nau::scene;

Pass *
PAssCudaPI::Create(const std::string &passName) {

	return new PAssCudaPI(passName);
}


PAssCudaPI::PAssCudaPI(const std::string &passName) :
	Pass (passName) {

	m_ClassName = "cudaPI";
}


PAssCudaPI::~PAssCudaPI(void) {

}

int first = 0;
int ranInit = 0;


struct cudaGraphicsResource* cuda_ssbo_Position;
float4* dptrssboPosition;

int* dptrssboIndex;
int* dptrssboTempIndex;
float4* dptrssboVelocity;
int* dptrssboCellStart;
int* dptrssboCellEnd;
float* dptrssboDensity;
float* dptrssboPressure;
int* dptrssboAdj;
float4* dptrssboForce;
float4* dptrssboGravity;
float4* dptrssboSurfaceNormal;
float4* dptrssboSurfaceTension;
float4* dptrssboViscosity;
float4* dptrssboAcceleration;


#include "mySort.h"
void
PAssCudaPI::prepare (void) {

	int densityPressure = 1;
	printf("Start frame -----------------");

	if (first == 0) {
		
		IBuffer* bPosition = RESOURCEMANAGER->getBuffer("simulationLib::Position");
		int buffIdPosition = bPosition->getPropi(IBuffer::ID);
		
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Position, buffIdPosition, cudaGraphicsMapFlagsNone);
		
		cudaMalloc((void**)&dptrssboIndex, 216000 * sizeof(int));
		cudaMalloc((void**)&dptrssboTempIndex, 216000 * sizeof(int));
		cudaMalloc((void **)&dptrssboVelocity, 216000 * sizeof(float4));
		cudaMalloc((void **)&dptrssboCellStart, 2000000 * sizeof(int));
		cudaMalloc((void **)&dptrssboCellEnd, 2000000 * sizeof(int));
		cudaMalloc((void **)&dptrssboDensity, 216000 * sizeof(float));
		cudaMalloc((void **)&dptrssboPressure, 216000 * sizeof(float));
		cudaMalloc((void **)&dptrssboAdj, 216000*500 * sizeof(int));
		cudaMalloc((void **)&dptrssboForce, 216000 * sizeof(float4));
		cudaMalloc((void **)&dptrssboGravity, 216000 * sizeof(float4));
		cudaMalloc((void **)&dptrssboSurfaceNormal, 216000 * sizeof(float4));
		cudaMalloc((void **)&dptrssboSurfaceTension, 216000 * sizeof(float4));
		cudaMalloc((void **)&dptrssboViscosity, 216000 * sizeof(float4));
		cudaMalloc((void **)&dptrssboAcceleration, 216000 * sizeof(float4));

		cudaMemset(dptrssboIndex,0, 216000 * sizeof(int));
		cudaMemset(dptrssboTempIndex, 0, 216000 * sizeof(int));
		cudaMemset(dptrssboVelocity, 0, 216000 * sizeof(float4));
		
		cudaMemset(dptrssboDensity, 0, 216000 * sizeof(float));
		cudaMemset(dptrssboPressure, 0, 216000 * sizeof(float));
		cudaMemset(dptrssboAdj, 0, 216000 * 500 * sizeof(int));
		cudaMemset(dptrssboForce, 0, 216000 * sizeof(float4));
		cudaMemset(dptrssboGravity, 0, 216000 * sizeof(float4));
		cudaMemset(dptrssboSurfaceNormal, 0, 216000 * sizeof(float4));
		cudaMemset(dptrssboSurfaceTension, 0, 216000 * sizeof(float4));
		cudaMemset(dptrssboViscosity, 0, 216000 * sizeof(float4));
		cudaMemset(dptrssboAcceleration, 0, 216000 * sizeof(float4));
		
		
		first = 1;
	}
	cudaDeviceSynchronize();
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	
	cudaMemset(dptrssboCellStart, 0, 2000000 * sizeof(int));
	cudaMemset(dptrssboCellEnd, 0, 2000000 * sizeof(int));
	
	cudaGraphicsMapResources(1, &cuda_ssbo_Position, NULL);

	cudaDeviceSynchronize();
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "cudaMemset = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

	
	size_t num_bytesssbo_Position;
	
	
	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboPosition, &num_bytesssbo_Position, cuda_ssbo_Position);
	
	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "cudaGraphicsResourceGetMappedPointer = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

	if (ranInit == 0) {
		
		cudaInit(dptrssboIndex, dptrssboTempIndex, dptrssboPosition);
		
		ranInit = 1;

	}
	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	
	mysort(dptrssboIndex, dptrssboPosition, dptrssboTempIndex, dptrssboVelocity, 216000);
	
	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "mysort = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
	//calls a kernel that counds each index to cellstart and cellend
	
	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	kernelWraper(dptrssboIndex, dptrssboCellStart, dptrssboCellEnd);

	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "Cell start/end = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;


	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();

	cudaDensityPressure(dptrssboPosition, dptrssboIndex, dptrssboCellStart, dptrssboCellEnd, dptrssboDensity, dptrssboPressure, dptrssboAdj);
	
	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "cudaDensityPressure = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
	
	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	
	cudaForce(dptrssboPosition, dptrssboVelocity, dptrssboForce, dptrssboDensity, dptrssboPressure, dptrssboGravity, dptrssboSurfaceNormal, dptrssboSurfaceTension, dptrssboViscosity, dptrssboAdj);
	
	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "cudaForce = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	cudaIntegrate(dptrssboPosition, dptrssboVelocity, dptrssboForce, dptrssboDensity, dptrssboGravity, dptrssboSurfaceTension, dptrssboViscosity, dptrssboAcceleration);
	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "Integrate = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	cudaUpdateIndex(dptrssboPosition, dptrssboIndex, dptrssboTempIndex);
	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "UpdateIndex = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
	
	
	//begin = std::chrono::steady_clock::now();
	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	cudaGraphicsUnmapResources(1, &cuda_ssbo_Position, NULL);
	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "unmap = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
	printf("End frame -----------------");
	cudaDeviceSynchronize();
}


void
PAssCudaPI::restore (void) {

}


void 
PAssCudaPI::doPass (void) {


}

