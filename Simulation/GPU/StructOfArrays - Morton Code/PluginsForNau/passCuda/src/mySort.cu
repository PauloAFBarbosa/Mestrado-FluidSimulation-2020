#include "mySort.h"
#include <cuda.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>




void mysort(int ** dptrssbo){
	
	thrust::device_ptr<int> buff = thrust::device_pointer_cast(*(dptrssbo));
	thrust::stable_sort(buff, buff + 20000000);
}

