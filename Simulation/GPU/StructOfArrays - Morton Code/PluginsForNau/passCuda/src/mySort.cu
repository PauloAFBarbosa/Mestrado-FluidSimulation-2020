#include "mySort.h"
#include <cuda.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>




void mysort(int ** index1, float4 ** values1, int** index2, float4 ** values2,int particles){
	
	thrust::device_ptr<int> i1buff = thrust::device_pointer_cast(*(index1));
	thrust::device_ptr<float4> v1buff = thrust::device_pointer_cast(*(values1));
	thrust::device_ptr<int> i2buff = thrust::device_pointer_cast(*(index2));
	thrust::device_ptr<float4> v2buff = thrust::device_pointer_cast(*(values2));


	thrust::stable_sort_by_key(i1buff, i1buff + particles,v1buff);
	//cudaThreadSynchronize();
	thrust::stable_sort_by_key(i2buff, i2buff + particles, v2buff);

	//cudaThreadSynchronize();
}

