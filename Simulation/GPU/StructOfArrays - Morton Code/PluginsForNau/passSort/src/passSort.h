#ifndef DEPTHMAPPASS_H
#define DEPTHMAPPASS_H


#include "iNau.h"
#include "nau/render/pass.h"

class PassSort : public Pass
		{
		protected:
			std::shared_ptr<nau::scene::Camera> m_LightCamera;

		public:

			static Pass *Create(const std::string &passName);
			PassSort(const std::string &name);
			~PassSort(void);

			virtual void prepare (void);
			virtual void doPass (void);
			virtual void restore (void);
			
};

extern "C" {
#ifdef WIN32
	__declspec(dllexport) void *createPass(const char *s);
	__declspec(dllexport) void init(void *inau);
	__declspec(dllexport) char *getClassName();
#else
	void *createPass(const char *s);
	void init(void *inau);
	char *getClassName();
#endif	
}

#endif //DEPTHMAPPASS_H
