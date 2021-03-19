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

struct cudaGraphicsResource* cuda_ssbo_Index;
struct cudaGraphicsResource* cuda_ssbo_TempIndex;
struct cudaGraphicsResource* cuda_ssbo_Position;
struct cudaGraphicsResource* cuda_ssbo_Velocity;
struct cudaGraphicsResource* cuda_ssbo_CellStart;
struct cudaGraphicsResource* cuda_ssbo_CellEnd;

struct cudaGraphicsResource* cuda_ssbo_Density;
struct cudaGraphicsResource* cuda_ssbo_Pressure;
struct cudaGraphicsResource* cuda_ssbo_Adj;
struct cudaGraphicsResource* cuda_ssbo_Force;
struct cudaGraphicsResource* cuda_ssbo_Gravity;
struct cudaGraphicsResource* cuda_ssbo_SurfaceNormal;
struct cudaGraphicsResource* cuda_ssbo_SurfaceTension;
struct cudaGraphicsResource* cuda_ssbo_Viscosity;


#include "mySort.h"
void
PAssCudaPI::prepare (void) {

	int densityPressure = 0;

	if (first == 0) {
		//index
		IBuffer* bIndex = RESOURCEMANAGER->getBuffer("simulationLib::Index");
		int buffIdIndex = bIndex->getPropi(IBuffer::ID);
		//Position
		IBuffer* bPosition = RESOURCEMANAGER->getBuffer("simulationLib::Position");
		int buffIdPosition = bPosition->getPropi(IBuffer::ID);
		//TempIndex
		IBuffer* bTempIndex = RESOURCEMANAGER->getBuffer("simulationLib::TempIndex");
		int buffIdTempIndex = bTempIndex->getPropi(IBuffer::ID);
		//Velocity
		IBuffer* bVelocity = RESOURCEMANAGER->getBuffer("simulationLib::Velocity");
		int buffIdVelocity = bVelocity->getPropi(IBuffer::ID);
		//cell start and end
		IBuffer* bCellStart = RESOURCEMANAGER->getBuffer("simulationLib::CellStart");
		int buffIdCellStart = bCellStart->getPropi(IBuffer::ID);
		IBuffer* bCellEnd = RESOURCEMANAGER->getBuffer("simulationLib::CellEnd");
		int buffIdCellEnd = bCellEnd->getPropi(IBuffer::ID);

		// register this buffer object with CUDA
		//So devia chamar isto uma vez
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Index, buffIdIndex, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_TempIndex, buffIdTempIndex, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Position, buffIdPosition, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Velocity, buffIdVelocity, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_CellStart, buffIdCellStart, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_CellEnd, buffIdCellEnd, cudaGraphicsMapFlagsNone);

		if (densityPressure == 1) {
			IBuffer* bDensity = RESOURCEMANAGER->getBuffer("simulationLib::Density");
			int buffIdDensity = bDensity->getPropi(IBuffer::ID);
			IBuffer* bPressure = RESOURCEMANAGER->getBuffer("simulationLib::Pressure");
			int buffIdPressure = bPressure->getPropi(IBuffer::ID);
			IBuffer* bAdj = RESOURCEMANAGER->getBuffer("simulationLib::Adj");
			int buffIdAdj = bAdj->getPropi(IBuffer::ID);

			IBuffer* bForce = RESOURCEMANAGER->getBuffer("simulationLib::Force");
			int buffIdForce = bForce->getPropi(IBuffer::ID);
			IBuffer* bGravity = RESOURCEMANAGER->getBuffer("simulationLib::Gravity");
			int buffIdGravity = bGravity->getPropi(IBuffer::ID);
			IBuffer* bSurfaceNormal = RESOURCEMANAGER->getBuffer("simulationLib::SurfaceNormal");
			int buffIdSurfaceNormal = bSurfaceNormal->getPropi(IBuffer::ID);
			IBuffer* bSurfaceTension = RESOURCEMANAGER->getBuffer("simulationLib::SurfaceTension");
			int buffIdSurfaceTension = bSurfaceTension->getPropi(IBuffer::ID);
			IBuffer* bViscosity = RESOURCEMANAGER->getBuffer("simulationLib::Viscosity");
			int buffIdViscosity = bViscosity->getPropi(IBuffer::ID);

			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Density, buffIdDensity, cudaGraphicsMapFlagsNone);
			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Pressure, buffIdPressure, cudaGraphicsMapFlagsNone);
			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Adj, buffIdAdj, cudaGraphicsMapFlagsNone);

			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Force, buffIdForce, cudaGraphicsMapFlagsNone);
			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Gravity, buffIdGravity, cudaGraphicsMapFlagsNone);
			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_SurfaceNormal, buffIdSurfaceNormal, cudaGraphicsMapFlagsNone);
			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_SurfaceTension, buffIdSurfaceTension, cudaGraphicsMapFlagsNone);
			cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Viscosity, buffIdViscosity, cudaGraphicsMapFlagsNone);
		}
		
		first = 1;
	}
	

	// map OpenGL buffer object for writing from CUDA
	int* dptrssboIndex;
	int* dptrssboTempIndex;
	float4 * dptrssboPosition;
	float4 * dptrssboVelocity;
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

	cudaDeviceSynchronize();
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	
	cudaGraphicsMapResources(1, &cuda_ssbo_Index, NULL);
	cudaGraphicsMapResources(1, &cuda_ssbo_TempIndex, NULL);
	cudaGraphicsMapResources(1, &cuda_ssbo_Position, NULL);
	cudaGraphicsMapResources(1, &cuda_ssbo_Velocity, NULL);
	cudaGraphicsMapResources(1, &cuda_ssbo_CellStart, NULL);
	cudaGraphicsMapResources(1, &cuda_ssbo_CellEnd, NULL);
	
	if (densityPressure == 1) {
		cudaGraphicsMapResources(1, &cuda_ssbo_Density, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_Pressure, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_Adj, NULL);

		cudaGraphicsMapResources(1, &cuda_ssbo_Force, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_Gravity, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_SurfaceNormal, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_SurfaceTension, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_Viscosity, NULL);
	}

	cudaDeviceSynchronize();
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "cudaGraphicsMapResources = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

	size_t num_bytesssbo_Index;
	size_t num_bytesssbo_TempIndex;
	size_t num_bytesssbo_Position;
	size_t num_bytesssbo_Velocity;
	size_t num_bytesssbo_CellStart;
	size_t num_bytesssbo_CellEnd;
	
	size_t num_bytesssbo_Density;
	size_t num_bytesssbo_Pressure;
	size_t num_bytesssbo_Adj;

	size_t num_bytesssbo_Force;
	size_t num_bytesssbo_Gravity;
	size_t num_bytesssbo_SurfaceNormal;
	size_t num_bytesssbo_SurfaceTension;
	size_t num_bytesssbo_Viscosity;
	
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboIndex, &num_bytesssbo_Index, cuda_ssbo_Index);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboTempIndex, &num_bytesssbo_TempIndex, cuda_ssbo_TempIndex);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboPosition, &num_bytesssbo_Position, cuda_ssbo_Position);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboVelocity, &num_bytesssbo_Velocity, cuda_ssbo_Velocity);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboCellStart, &num_bytesssbo_CellStart, cuda_ssbo_CellStart);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboCellEnd, &num_bytesssbo_CellEnd, cuda_ssbo_CellEnd);

	if (densityPressure == 1) {
		
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboDensity, &num_bytesssbo_Density, cuda_ssbo_Density);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboPressure, &num_bytesssbo_Pressure, cuda_ssbo_Pressure);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboAdj, &num_bytesssbo_Adj, cuda_ssbo_Adj);

		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboForce, &num_bytesssbo_Force, cuda_ssbo_Force);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboGravity, &num_bytesssbo_Gravity, cuda_ssbo_Gravity);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboSurfaceNormal, &num_bytesssbo_SurfaceNormal, cuda_ssbo_SurfaceNormal);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboSurfaceTension, &num_bytesssbo_SurfaceTension, cuda_ssbo_SurfaceTension);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboViscosity, &num_bytesssbo_Viscosity, cuda_ssbo_Viscosity);
	}

	mysort(dptrssboIndex,dptrssboPosition, dptrssboTempIndex, dptrssboVelocity,216000);

	//calls a kernel that counds each index to cellstart and cellend
	kernelWraper(dptrssboIndex, dptrssboCellStart, dptrssboCellEnd);

	if (densityPressure == 1) {
		cudaDeviceSynchronize();
		begin = std::chrono::steady_clock::now();
		cudaDensityPressure(dptrssboPosition, dptrssboIndex, dptrssboCellStart, dptrssboCellEnd, dptrssboDensity, dptrssboPressure, dptrssboAdj);
		cudaDeviceSynchronize();
		end = std::chrono::steady_clock::now();
		std::cout << "DensityPressure = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

		cudaForce(dptrssboPosition, dptrssboVelocity, dptrssboForce, dptrssboDensity, dptrssboPressure, dptrssboGravity, dptrssboSurfaceNormal, dptrssboSurfaceTension, dptrssboViscosity, dptrssboAdj);
	}
	
	cudaDeviceSynchronize();
	begin = std::chrono::steady_clock::now();
	cudaGraphicsUnmapResources(1, &cuda_ssbo_Index, NULL);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_TempIndex, NULL);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_Position, NULL);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_Velocity, NULL);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_CellStart, NULL);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_CellEnd, NULL);

	if (densityPressure == 1) {
		cudaGraphicsUnmapResources(1, &cuda_ssbo_Density, NULL);
		cudaGraphicsUnmapResources(1, &cuda_ssbo_Pressure, NULL);
		cudaGraphicsUnmapResources(1, &cuda_ssbo_Adj, NULL);
	}

	cudaDeviceSynchronize();
	end = std::chrono::steady_clock::now();
	std::cout << "cudaGraphicsUnmapResources = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
	
}


void
PAssCudaPI::restore (void) {

}


void 
PAssCudaPI::doPass (void) {


}

