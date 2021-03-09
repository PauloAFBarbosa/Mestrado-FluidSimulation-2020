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

#include "mySort.h"
void
PAssCudaPI::prepare (void) {

	//vou buscar o identificador do buffer
	IBuffer* b = RESOURCEMANAGER->getBuffer("particleLib::Index");

	int buffId = b->getPropi(IBuffer::ID);

	struct cudaGraphicsResource* cuda_ssbo_resource;

	void* d_ssbo_buffer = NULL;

	// register this buffer object with CUDA
	cudaGraphicsGLRegisterBuffer(&cuda_ssbo_resource, buffId, cudaGraphicsMapFlagsNone);

	// map OpenGL buffer object for writing from CUDA
	int* dptrssbo;

	cudaGraphicsMapResources(1, &cuda_ssbo_resource, 0);

	size_t num_bytesssbo;

	cudaGraphicsResourceGetMappedPointer((void**)&dptrssbo, &num_bytesssbo, cuda_ssbo_resource);

	mysort(&dptrssbo);

	cudaGraphicsUnmapResources(0, &cuda_ssbo_resource, 0);

}


void
PAssCudaPI::restore (void) {

}


void 
PAssCudaPI::doPass (void) {


}

