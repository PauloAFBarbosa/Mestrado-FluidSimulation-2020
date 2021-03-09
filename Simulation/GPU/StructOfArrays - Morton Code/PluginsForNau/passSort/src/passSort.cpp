#include "passSort.h"

#include "iNau.h"
#include "nau.h"
#include "nau/geometry/frustum.h"
#include "nau/render/passFactory.h"
#include "nau/material/iTexture.h"
#include "../../build/plugins/passSort/mySort.cu"

#include <glbinding/gl/gl.h>
#include <glbinding/Binding.h>

#ifdef _WIN32
#  define WINDOWS_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#endif


#include <cuda_gl_interop.h>





//static nau::INau *NAU_INTERFACE;

static char className[] = "passSort";

#ifdef WIN32
#define EXPORT __declspec(dllexport)
#elif __linux__
#define EXPORT extern "C"
#endif

EXPORT
void*
createPass(const char* s) {

	std::shared_ptr<PassSort>* p = new std::shared_ptr<PassSort>(new PassSort(s));
	return p;
}


EXPORT
void
init(void* nauInst) {

	//NAU_INTERFACE = (nau::INau *)nauInst;
	INau::SetInterface((nau::INau*)nauInst);
	nau::Nau::SetInstance((nau::Nau*)nauInst);
	glbinding::Binding::initialize(false);
}


EXPORT
char*
getClassName() {

	return className;
}


using namespace nau::geometry;
using namespace nau::math;
using namespace nau::render;
using namespace nau::scene;

Pass*
PassSort::Create(const std::string& passName) {

	return new PassSort(passName);
}


PassSort::PassSort(const std::string& passName) :
	Pass(passName) {

	//m_ClassName = "passSort";
	//std::string camName = passName + "-LightCam";
	//m_Viewport[0] = RENDERMANAGER->createViewport(camName);
	//m_LightCamera = RENDERMANAGER->getCamera(camName);
}


PassSort::~PassSort(void) {

}

void
PassSort::prepare(void) {
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

	cudaGraphicsResourceGetMappedPointer((void**)&dptrssbo, &num_bytesssbo,cuda_ssbo_resource);

	mySort(dptrssbo);

	cudaGraphicsUnmapResources(0, &cuda_ssbo_resource, 0);
}


void
PassSort::restore(void) {

	
}


void
PassSort::doPass(void) {

}

