#include "passPICudaHash.h"

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

#include "nau/debug/profile.h"
//using namespace gl;

//static nau::INau *NAU_INTERFACE;

static char className[] = "cudaPIHash";

#ifdef WIN32
#define EXPORT __declspec(dllexport)
#elif __linux__
#define EXPORT extern "C"
#endif

EXPORT 
void *
createPass(const char *s) {

	std::shared_ptr<PAssCudaPIHash> *p = new std::shared_ptr<PAssCudaPIHash>(new PAssCudaPIHash(s));
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
PAssCudaPIHash::Create(const std::string &passName) {

	return new PAssCudaPIHash(passName);
}


PAssCudaPIHash::PAssCudaPIHash(const std::string &passName) :
	Pass (passName) {

	m_ClassName = "cudaPIHash";
}


PAssCudaPIHash::~PAssCudaPIHash(void) {

}

int first = 0;

struct cudaGraphicsResource* cuda_ssbo_Index;
struct cudaGraphicsResource* cuda_ssbo_TempIndex;
struct cudaGraphicsResource* cuda_ssbo_Position;
struct cudaGraphicsResource* cuda_ssbo_Velocity;
struct cudaGraphicsResource* cuda_ssbo_CellStart;
struct cudaGraphicsResource* cuda_ssbo_CellEnd;

int particleNumber=0;
int hashSize = 0;

#include "mySortHash.h"
void
PAssCudaPIHash::prepare (void) {

	//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	if (first == 0) {
		

		//ir buscar o nr de particulas aqui
		void* d = (Data*)NAU->getAttributeValue("RENDERER", "CURRENT", "Number_Particles", 0);
		particleNumber = *(int*)d;

		//ir buscar o nr de particulas aqui
		d = (Data*)NAU->getAttributeValue("RENDERER", "CURRENT", "HASHSIZE", 0);
		hashSize = *(int*)d;

		//index
		IBuffer* bIndex = RESOURCEMANAGER->getBuffer("simulationLib::IndexPosition");
		int buffIdIndex = bIndex->getPropi(IBuffer::ID);
		//Position
		IBuffer* bPosition = RESOURCEMANAGER->getBuffer("simulationLib::Position");
		int buffIdPosition = bPosition->getPropi(IBuffer::ID);
		//TempIndex
		IBuffer* bTempIndex = RESOURCEMANAGER->getBuffer("simulationLib::IndexVelocity");
		int buffIdTempIndex = bTempIndex->getPropi(IBuffer::ID);
		//Velocity
		IBuffer* bVelocity = RESOURCEMANAGER->getBuffer("simulationLib::Velocity");
		int buffIdVelocity = bVelocity->getPropi(IBuffer::ID);
		//cell start and end

		IBuffer* bCellStart = RESOURCEMANAGER->getBuffer("simulationLib::Offsets");
		int buffIdCellStart = bCellStart->getPropi(IBuffer::ID);
		IBuffer* bCellEnd = RESOURCEMANAGER->getBuffer("simulationLib::BucketSizes");
		int buffIdCellEnd = bCellEnd->getPropi(IBuffer::ID);

		// register this buffer object with CUDA
		//So devia chamar isto uma vez
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Index, buffIdIndex, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_TempIndex, buffIdTempIndex, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Position, buffIdPosition, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Velocity, buffIdVelocity, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_CellStart, buffIdCellStart, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_CellEnd, buffIdCellEnd, cudaGraphicsMapFlagsNone);
		
		
		first = 1;
	}
	

	// map OpenGL buffer object for writing from CUDA
	int* dptrssboIndex;
	int* dptrssboTempIndex;
	float4 * dptrssboPosition;
	float4 * dptrssboVelocity;
	int* dptrssboCellStart;
	int* dptrssboCellEnd;

	//cudaDeviceSynchronize();
	//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	{
		PROFILE("MapResources");
		cudaGraphicsMapResources(1, &cuda_ssbo_Index, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_TempIndex, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_Position, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_Velocity, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_CellStart, NULL);
		cudaGraphicsMapResources(1, &cuda_ssbo_CellEnd, NULL);
	}

	//cudaDeviceSynchronize();
	//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	//std::cout << "cudaGraphicsMapResources = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

	size_t num_bytesssbo_Index;
	size_t num_bytesssbo_TempIndex;
	size_t num_bytesssbo_Position;
	size_t num_bytesssbo_Velocity;
	size_t num_bytesssbo_CellStart;
	size_t num_bytesssbo_CellEnd;
	
	{
		PROFILE("GetMappedPointer");
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboIndex, &num_bytesssbo_Index, cuda_ssbo_Index);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboTempIndex, &num_bytesssbo_TempIndex, cuda_ssbo_TempIndex);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboPosition, &num_bytesssbo_Position, cuda_ssbo_Position);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboVelocity, &num_bytesssbo_Velocity, cuda_ssbo_Velocity);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboCellStart, &num_bytesssbo_CellStart, cuda_ssbo_CellStart);
		cudaGraphicsResourceGetMappedPointer((void**)&dptrssboCellEnd, &num_bytesssbo_CellEnd, cuda_ssbo_CellEnd);
	}
	
	{
		PROFILE("MySortThrust");
		mysort(dptrssboIndex, dptrssboPosition, dptrssboTempIndex, dptrssboVelocity, particleNumber);
	}
	
	{
		PROFILE("CellStartEnd");
		//calls a kernel that counds each index to cellstart and cellend
		kernelWraper(dptrssboIndex, dptrssboCellStart, dptrssboCellEnd, particleNumber, hashSize);
	}
	
	
	//cudaDeviceSynchronize();
	//begin = std::chrono::steady_clock::now();
	{
		PROFILE("UnmapResources");
		cudaGraphicsUnmapResources(1, &cuda_ssbo_Index, NULL);
		cudaGraphicsUnmapResources(1, &cuda_ssbo_TempIndex, NULL);
		cudaGraphicsUnmapResources(1, &cuda_ssbo_Position, NULL);
		cudaGraphicsUnmapResources(1, &cuda_ssbo_Velocity, NULL);
		cudaGraphicsUnmapResources(1, &cuda_ssbo_CellStart, NULL);
		cudaGraphicsUnmapResources(1, &cuda_ssbo_CellEnd, NULL);
	}
	//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	
	//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[milliseconds]" << std::endl;
	
}


void
PAssCudaPIHash::restore (void) {

}


void 
PAssCudaPIHash::doPass (void) {


}

