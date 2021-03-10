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

#include "mySort.h"
void
PAssCudaPI::prepare (void) {

	//vou buscar o identificador do buffer
	//preciso do Index,Position
	//preciso do TempIndex,Velocity
	if (first == 0) {
		//index
		IBuffer* bIndex = RESOURCEMANAGER->getBuffer("particleLib::Index");
		int buffIdIndex = bIndex->getPropi(IBuffer::ID);
		//Position
		IBuffer* bPosition = RESOURCEMANAGER->getBuffer("particleLib::Position");
		int buffIdPosition = bPosition->getPropi(IBuffer::ID);
		//TempIndex
		IBuffer* bTempIndex = RESOURCEMANAGER->getBuffer("particleLib::TempIndex");
		int buffIdTempIndex = bTempIndex->getPropi(IBuffer::ID);
		//Velocity
		IBuffer* bVelocity = RESOURCEMANAGER->getBuffer("particleLib::Velocity");
		int buffIdVelocity = bVelocity->getPropi(IBuffer::ID);

		// register this buffer object with CUDA
		//So devia chamar isto uma vez
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Index, buffIdIndex, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_TempIndex, buffIdTempIndex, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Position, buffIdPosition, cudaGraphicsMapFlagsNone);
		cudaGraphicsGLRegisterBuffer(&cuda_ssbo_Velocity, buffIdVelocity, cudaGraphicsMapFlagsNone);
		first = 1;
	}
	

	// map OpenGL buffer object for writing from CUDA
	int* dptrssboIndex;
	int* dptrssboTempIndex;
	float4 * dptrssboPosition;
	float4 * dptrssboVelocity;

	cudaGraphicsMapResources(1, &cuda_ssbo_Index, 0);
	cudaGraphicsMapResources(1, &cuda_ssbo_TempIndex, 0);
	cudaGraphicsMapResources(1, &cuda_ssbo_Position, 0);
	cudaGraphicsMapResources(1, &cuda_ssbo_Velocity, 0);

	size_t num_bytesssbo_Index;
	size_t num_bytesssbo_TempIndex;
	size_t num_bytesssbo_Position;
	size_t num_bytesssbo_Velocity;

	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboIndex, &num_bytesssbo_Index, cuda_ssbo_Index);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboTempIndex, &num_bytesssbo_TempIndex, cuda_ssbo_TempIndex);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboPosition, &num_bytesssbo_Position, cuda_ssbo_Position);
	cudaGraphicsResourceGetMappedPointer((void**)&dptrssboVelocity, &num_bytesssbo_Velocity, cuda_ssbo_Velocity);

	mysort(&dptrssboIndex,&dptrssboPosition, &dptrssboTempIndex, &dptrssboVelocity,216000);

	cudaGraphicsUnmapResources(1, &cuda_ssbo_Index, 0);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_TempIndex, 0);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_Position, 0);
	cudaGraphicsUnmapResources(1, &cuda_ssbo_Velocity, 0);

}


void
PAssCudaPI::restore (void) {

}


void 
PAssCudaPI::doPass (void) {


}

