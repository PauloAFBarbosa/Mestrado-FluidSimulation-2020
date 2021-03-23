
#include <cuda.h>

void mysort(int* index1, float4 * values1, int* index2, float4 * values2,int particles);

void kernelWraper(int* dptrssboIndex, int* dptrssboCellStart, int* dptrssboCellEnd);

void cudaDensityPressure(float4* position,int* dptrssboIndex, int* dptrssboCellStart, int* dptrssboCellEnd, float* dptrssboDensity, float* dptrssboPressure, int * dptrssboAdj);

void cudaForce(float4 * dptrssboPosition, float4 * dptrssboVelocity, float4 * dptrssboForce, float * dptrssboDensity, float * dptrssboPressure, float4 * dptrssboGravity, float4 * dptrssboSurfaceNormal, float4 * dptrssboSurfaceTension, float4 * dptrssboViscosity, int * dptrssboAdj);

void cudaIntegrate(float4* dptrssboPosition, float4* dptrssboVelocity, float4* dptrssboForce, float* dptrssboDensity, float4* dptrssboGravity, float4* dptrssboSurfaceTension, float4* dptrssboViscosity, float4* dptrssboAcceleration);

void cudaUpdateIndex(float4* dptrssboPosition, int* dptrssboIndex, int* dptrssboTempIndex);

void cudaInit(int* dptrssboIndex, int* dptrssboTempIndex, float4* dptrssboPosition);