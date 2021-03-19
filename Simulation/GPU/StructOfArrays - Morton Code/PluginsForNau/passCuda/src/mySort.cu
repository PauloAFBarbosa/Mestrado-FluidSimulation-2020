#include "mySort.h"
#include <cuda.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <chrono>
#include <cuda_runtime.h>
#include "cutil_math.h"
//#include <math_functions.h>



void mysort(int * index1, float4 * values1, int* index2, float4 * values2,int particles){
	
	
	
	thrust::device_ptr<int> i1buff = thrust::device_pointer_cast((index1));
	thrust::device_ptr<float4> v1buff = thrust::device_pointer_cast((values1));
	thrust::device_ptr<int> i2buff = thrust::device_pointer_cast((index2));
	thrust::device_ptr<float4> v2buff = thrust::device_pointer_cast((values2));

	
	thrust::sort_by_key(i1buff, i1buff + particles, v1buff);
	thrust::sort_by_key(i2buff, i2buff + particles, v2buff);
	
}

//functions for density Pressure
//Morton code --------------------------------------------


__device__
unsigned int part1by2(unsigned int n) {
    n &= 0x000003ff;
    n = (n ^ (n << 16)) & 0xff0000ff;
    n = (n ^ (n << 8)) & 0x0300f00f;
    n = (n ^ (n << 4)) & 0x030c30c3;
    n = (n ^ (n << 2)) & 0x09249249;
    return n;
}

__device__
unsigned int unpart1by2(unsigned int n) {
    n &= 0x09249249;
    n = (n ^ (n >> 2)) & 0x030c30c3;
    n = (n ^ (n >> 4)) & 0x0300f00f;
    n = (n ^ (n >> 8)) & 0xff0000ff;
    n = (n ^ (n >> 16)) & 0x000003ff;
    return n;
}

__device__
unsigned int interleave3(unsigned int x, unsigned int y, unsigned int z) {
    return part1by2(x) | (part1by2(y) << 1) | (part1by2(z) << 2);
}

__device__
void deinterleave3(unsigned int n, unsigned int x, unsigned int y, unsigned int z) {
    x = unpart1by2(n);
    y = unpart1by2(n >> 1);
    z = unpart1by2(n >> 2);
}


__device__
bool contains(unsigned int arr[27], int size, unsigned int member) {
    bool ret = false;

    for (int i = 0; i < size; i++)
    {
        if (arr[i] == member)
            ret = true;
    }
    return ret;
}
__device__
int getAdj(float4 pos, float H,unsigned int ret[27]) {

    int retSize = 0;
    unsigned int offset = 43;
    unsigned int morton_x = unsigned int((pos.x / H) + offset);
    unsigned int morton_y = unsigned int((pos.y / H) + offset);
    unsigned int morton_z = unsigned int((pos.z / H) + offset);

    unsigned int morton_cell;

    if (morton_x > 0 && morton_y > 0 && morton_z > 0) {
        morton_cell = interleave3(morton_x - 1, morton_y - 1, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_y > 0 && morton_z > 0) {
        morton_cell = interleave3(morton_x, morton_y - 1, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_x > 0 && morton_z > 0) {
        morton_cell = interleave3(morton_x - 1, morton_y, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_z > 0) {
        morton_cell = interleave3(morton_x, morton_y, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_x > 0 && morton_y > 0) {
        morton_cell = interleave3(morton_x - 1, morton_y - 1, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_y > 0) {
        morton_cell = interleave3(morton_x, morton_y - 1, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_x > 0) {
        morton_cell = interleave3(morton_x - 1, morton_y, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }

    morton_cell = interleave3(morton_x, morton_y, morton_z);
    ret[retSize] = morton_cell;
    retSize++;

    //1864184 é o numero maximo que o morton code pode devolver num cubo de -2 a 2 
    if (morton_x < 1864184 && morton_y > 0 && morton_z > 0) {
        morton_cell = interleave3(morton_x + 1, morton_y - 1, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_x < 1864184 && morton_z > 0) {
        morton_cell = interleave3(morton_x + 1, morton_y, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_x < 1864184 && morton_y > 0) {
        morton_cell = interleave3(morton_x + 1, morton_y - 1, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x < 1864184) {
        morton_cell = interleave3(morton_x + 1, morton_y, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x > 0 && morton_y < 1864184 && morton_z > 0) {
        morton_cell = interleave3(morton_x - 1, morton_y + 1, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_y < 1864184 && morton_z > 0) {
        morton_cell = interleave3(morton_x, morton_y + 1, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_y < 1864184 && morton_x > 0) {
        morton_cell = interleave3(morton_x - 1, morton_y + 1, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_y < 1864184) {
        morton_cell = interleave3(morton_x, morton_y + 1, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x < 1864184 && morton_y < 1864184 && morton_z > 0) {
        morton_cell = interleave3(morton_x + 1, morton_y + 1, morton_z - 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x < 1864184 && morton_y < 1864184) {
        morton_cell = interleave3(morton_x + 1, morton_y + 1, morton_z);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x > 0 && morton_y > 0 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x - 1, morton_y - 1, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_y > 0 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x, morton_y - 1, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_x > 0 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x - 1, morton_y, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_z < 1864184) {
        morton_cell = interleave3(morton_x, morton_y, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x < 1864184 && morton_y >0 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x + 1, morton_y - 1, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }
    if (morton_x < 1864184 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x + 1, morton_y, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x > 0 && morton_y < 1864184 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x - 1, morton_y + 1, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_y < 1864184 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x, morton_y + 1, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    if (morton_x < 1864184 && morton_y < 1864184 && morton_z < 1864184) {
        morton_cell = interleave3(morton_x + 1, morton_y + 1, morton_z + 1);
        ret[retSize] = morton_cell;
        retSize++;
    }

    return retSize;

}

__device__
float useDefaultKernel(float4 distVector, float supportRadius) {
    
    float dist = length(distVector);
    
    if (dist > supportRadius) {

        return 0.0;
    }
    else {
        //printf("Vizinho e vai devolver -> %f \n", (315 / (64 * 3.141592653589793 * pow(supportRadius, 9.0f))) * pow(supportRadius * supportRadius - dist * dist, 3.0f));
        return (315 / (64 * 3.141592653589793 * pow(supportRadius, 9.0f))) * pow(supportRadius * supportRadius - dist * dist, 3.0f);
    }
}
__device__
unsigned int hashFunction(float4 pos, double H, int size) {
    
    int p1 = 2693;
    int p2 = 3163;
    int p3 = 4091;

    int part1 = (int((pos.x / H)) * p1);
    int part2 = (int((pos.y / H)) * p2);
    int part3 = (int((pos.z / H)) * p3);

    unsigned int ret = unsigned int((part1 ^ part2 ^ part3) % size);

    return ret;
}
//----------------


__global__
void count(int * indexes , int* CellStart, int* CellEnd)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	int indexCount = indexes[i];
	atomicAdd(CellStart + indexCount, 1);

	atomicAdd(CellEnd + indexCount, 1);
	
}

void kernelWraper(int * dptrssboIndex, int* dptrssboCellStart, int* dptrssboCellEnd) {

	count<<<6750, 32 >>>(dptrssboIndex, dptrssboCellStart, dptrssboCellEnd);
	
	

	thrust::device_ptr<int> cellstartThrust = thrust::device_pointer_cast((dptrssboCellStart));
	thrust::exclusive_scan(cellstartThrust, cellstartThrust + 2000000, cellstartThrust);
	
}



__global__
void densityPressureKernel(float4 * dptrssboPosition,int* dptrssboIndex, int* dptrssboCellStart, int* dptrssboCellEnd, float* dptrssboDensity, float* dptrssboPressure, int* dptrssboAdj)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
    int conta = 0;

    unsigned int ret[27];
    int retSize = 0;
    float H = 0.0457;
    unsigned int offset = 43;
    unsigned int morton_x = unsigned int((dptrssboPosition[index].x / H) + offset);
    unsigned int morton_y = unsigned int((dptrssboPosition[index].y / H) + offset);
    unsigned int morton_z = unsigned int((dptrssboPosition[index].z / H) + offset);
    retSize = getAdj(dptrssboPosition[index], H, ret);

    //printf("Chegou aqui 1\n");
    
    // compute density
    float sum = 1;

    //mudar isto para um var local
    dptrssboAdj[index * 500] = 0;
    //printf("Chegou aqui 2\n");
    int vizinhos = 0;

    int naovizinhos = 0;

    for (int j = 0; j < retSize; j++)
    {
        //uint bucket = uint(adjMat[(28*b)+(j+1)]) ;
        unsigned int cell = ret[j];

        int from = dptrssboCellStart[cell];
        int to = from + dptrssboCellEnd[cell];
        //printf("Chegou aqui 3\n");
        for (int i = from; i < to; i++)
        {


            float4 arg;
            arg = dptrssboPosition[index] - dptrssboPosition[i];

            float ker_res = useDefaultKernel(arg, H);

            sum += 0.02 * ker_res;
            //printf("Chegou aqui 4\n");

            if (ker_res != 0) {

                conta++;
                int vizinhos = dptrssboAdj[index * 500];
                //Aqui é MAXADJ porque no maximo vai dar para guardar 99 vizinhos mais 1 posição para dizer quantos vizinhos se guardou
                dptrssboAdj[index * 500 + vizinhos + 1] = i;
                dptrssboAdj[index * 500]++;
                vizinhos++;
                //printf("Chegou aqui 5\n");
            }
            else
                naovizinhos++;


        }
    }

    //printf("Chegou aqui 6 vizinhos %d naovizinhos %d\n",vizinhos,naovizinhos);
    dptrssboDensity[index] = sum;

    //printf("Chegou aqui 7\n");
    // compute pressure
    dptrssboPressure[index] = 3.0 * (sum - 998.29);
}

__device__
float useDefaultKernel_laplacian(float4 distVector, float supportRadius) {
    float dist = length(distVector);
    if (dist > supportRadius)
        return 0.0f;
    else
        return -(945 / (32 * 3.141592653589793 * pow(supportRadius, 9.0f))) * (supportRadius * supportRadius - dist * dist) * (3 * supportRadius * supportRadius - 7 * dist * dist);
}
__device__
float useViscosityKernel_laplacian(float4 distVector, float supportRadius) {
    float dist = length(distVector);
    if (dist > supportRadius)
        return 0.0f;
    else
        return (45 / (3.141592653589793 * pow(supportRadius, 6.0f))) * (supportRadius - dist);
}
__device__
float4 useDefaultKernel_gradient(float4 distVector, float supportRadius) {
    float dist = length(distVector);
    if (dist > supportRadius) {
        return make_float4(0);
    }
    else {
        return -(distVector * (945 / (32 * 3.141592653589793 * pow(supportRadius, 9.0f))) * pow(supportRadius * supportRadius - dist * dist, 2.0f));
    }
}



__device__
float4 usePressureKernel_gradient(float4 distVector, float supportRadius) {
    float dist = length(distVector);
    if (dist > supportRadius) {

        return make_float4(0);

    }
    else
    {

        float4 normalized;
        normalized = normalize(distVector);

        return -(normalized * (45 / (3.141592653589793 * pow(supportRadius, 6.0f))) * pow(supportRadius - dist, 2.0f));

    }
}

__global__
void ForceKernel(float4* dptrssboPosition, float4* dptrssboVelocity, float4* dptrssboForce, float* dptrssboDensity, float* dptrssboPressure, float4* dptrssboGravity, float4* dptrssboSurfaceNormal, float4* dptrssboSurfaceTension, float4* dptrssboViscosity, int* dptrssboAdj)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int ret[27];
    int retSize = 0;
    float H = 0.0457;

    retSize = getAdj(dptrssboPosition[index], H, ret);

    float4 sumViscosity = make_float4(0);

    float4 sumForce= make_float4(0);
    

    float4 sumSurfaceNormal= make_float4(0);
    
    float sum = 0;
    float MASS = 0.02;

    int count;
    count = 0;

    int vizinhos = dptrssboAdj[index * 500];

    for (int i = 0; i < vizinhos; i++) {
        //Viscosity

        if (dptrssboPosition[index].x == dptrssboPosition[dptrssboAdj[index * 500 + i + 1]].x && dptrssboPosition[index].y == dptrssboPosition[dptrssboAdj[index * 500 + i + 1]].y && dptrssboPosition[index].z == dptrssboPosition[dptrssboAdj[index * 500 + i + 1]].z)
            continue;
        float4 arg;
        arg = dptrssboPosition[index] - dptrssboPosition[dptrssboAdj[index * 500 + i + 1]];

        sumViscosity += useViscosityKernel_laplacian(arg, H) * (dptrssboVelocity[dptrssboAdj[index * 500 + i + 1]] - dptrssboVelocity[index]) * (MASS / dptrssboDensity[dptrssboAdj[index * 500 + i + 1]]);

        //Force
  
        arg = dptrssboPosition[index] - dptrssboPosition[dptrssboAdj[index * 500 + i + 1]];

        
        //printf("ret Kernel %f %f %f \n", usePressureKernel_gradient(arg, H).x, usePressureKernel_gradient(arg, H).y, usePressureKernel_gradient(arg, H).z);
        sumForce += usePressureKernel_gradient(arg, H) * (dptrssboPressure[index] / (dptrssboDensity[index] * dptrssboDensity[index]) + dptrssboPressure[dptrssboAdj[index * 500 + i + 1]] / (dptrssboDensity[dptrssboAdj[index * 500 + i + 1]] * dptrssboDensity[dptrssboAdj[index * 500 + i + 1]])) * MASS;

        //SurfaceNormal
        arg = dptrssboPosition[index] - dptrssboPosition[dptrssboAdj[index * 500 + i + 1]];
        
        sumSurfaceNormal += useDefaultKernel_gradient(arg, H) * (MASS / dptrssboDensity[dptrssboAdj[index * 500 + i + 1]]);
    }

    dptrssboViscosity[index] = sumViscosity * 3.5;

    dptrssboForce[index] = -(sumForce * dptrssboDensity[index]);
    //debug[index] = vec4(density[index]);

    //tempPosition[index]=vec4(sumForce[0],sumForce[1],sumForce[2 ],0);

    dptrssboSurfaceNormal[index] = sumSurfaceNormal;

    dptrssboGravity[index] = make_float4(0.0,-9.8,0.0,0.0) * dptrssboDensity[index];



    if (length(dptrssboSurfaceNormal[index]) >= 7.065) {

        for (int i = 0; i < vizinhos; i++) {
            if (dptrssboPosition[index].x == dptrssboPosition[dptrssboAdj[index * 500 + i + 1]].x && dptrssboPosition[index].y == dptrssboPosition[dptrssboAdj[index * 500 + i + 1]].y && dptrssboPosition[index].z == dptrssboPosition[dptrssboAdj[index * 500 + i + 1]].z)
                continue;
            float4 arg;
            arg = dptrssboPosition[index] - dptrssboPosition[dptrssboAdj[index * 500 + i + 1]];
            
            
            sum += (MASS / dptrssboDensity[dptrssboAdj[index * 500 + i + 1]]) * useDefaultKernel_laplacian(arg, H);
        }




        float4 surfaceNormalNormalized;

        surfaceNormalNormalized = normalize(dptrssboSurfaceNormal[index]);

        dptrssboSurfaceTension[index] = -(surfaceNormalNormalized * 0.0728 * sum);

    }
    else {
        dptrssboSurfaceTension[index] =make_float4(0);
    }
}

void cudaDensityPressure(float4 * dptrssboPosition,int * dptrssboIndex,int * dptrssboCellStart,int * dptrssboCellEnd,float * dptrssboDensity, float * dptrssboPressure, int * dptrssboAdj) {
	densityPressureKernel << <6750, 32 >> > (dptrssboPosition,dptrssboIndex, dptrssboCellStart, dptrssboCellEnd, dptrssboDensity, dptrssboPressure, dptrssboAdj);
}

void cudaForce(float4* dptrssboPosition, float4* dptrssboVelocity, float4* dptrssboForce, float* dptrssboDensity, float* dptrssboPressure, float4* dptrssboGravity, float4* dptrssboSurfaceNormal, float4* dptrssboSurfaceTension, float4* dptrssboViscosity, int* dptrssboAdj) {
    ForceKernel << <6750, 32 >> > (dptrssboPosition, dptrssboVelocity, dptrssboForce, dptrssboDensity, dptrssboPressure, dptrssboGravity, dptrssboSurfaceNormal, dptrssboSurfaceTension, dptrssboViscosity, dptrssboAdj);
}
