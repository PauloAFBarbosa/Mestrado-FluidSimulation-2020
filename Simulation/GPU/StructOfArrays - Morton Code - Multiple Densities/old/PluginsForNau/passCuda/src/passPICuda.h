#ifndef CUDAPIPASS_H
#define CUDAPIPASS_H


#include "iNau.h"
#include "nau/render/pass.h"

class PAssCudaPI : public Pass
		{
		protected:

		public:

			static Pass *Create(const std::string &passName);
			PAssCudaPI(const std::string &name);
			~PAssCudaPI(void);

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

#endif //CUDAPIPASS_H
