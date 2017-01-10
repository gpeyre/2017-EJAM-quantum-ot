/*
 * (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
 * University, University College London. This code is released under the 
 * Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
 * license (http://creativecommons.org/licenses/by-nc-sa/4.0/).
 */

// nvcc --use_fast_math -O3 -ptx -arch=sm_30 cuda_feature.cu

__global__ void cuda_feature(unsigned int *perm, unsigned int *mins_idx, float *mins, float *fdists, unsigned int *A_feat, unsigned int *B_feat, float *A_col, float *B_col, int M2, int N, float w_col, int t)
{
/*
    int bIdx = gridDim.x*gridDim.y*blockIdx.z + gridDim.x*blockIdx.y + blockIdx.x;
    int idx = bIdx * blockDim.x * blockDim.y * blockDim.z 
            + blockDim.x * threadIdx.y
            + threadIdx.x;
*/

    //int idx = blockDim.x*gridDim.x*blockIdx.y + blockIdx.x * blockDim.x + threadIdx.x;
    int idx = blockDim.x * blockIdx.x + threadIdx.x;


    //fdists[idx] = 0;

    float dist = 0;

    unsigned int hamming = 0;
    for (int n = 0; n < N; n++)
	    hamming += __popc(A_feat[t + n*M2] ^ B_feat[idx + n*M2]);
    
    dist = (float)hamming / (float)(N*32);

    for (int n = 0; n < 3; n++)
	    dist += w_col*abs(A_col[t + n*M2] - B_col[idx + n*M2]);

    fdists[idx] = dist;

    __syncthreads();
    if (threadIdx.x != 0)
	    return;

    // Ugly recycling!
    dist = 999999;  // "infinity", eh.
    hamming = 0;

    for (int i = blockDim.x * blockIdx.x;
	 i < blockDim.x * (blockIdx.x+1); 
	 i++)
    {
	if (fdists[i] < dist)
	{
		dist = fdists[i];
		hamming = i;
	}
    }

    mins[blockIdx.x] = dist;
    mins_idx[blockIdx.x] = hamming;

    int prev = atomicInc(&mins_idx[gridDim.x], gridDim.x-1);
    if (idx != 0)
	    return;

    // Spin lock: wait until all blocks have incremented the
    // counter and it has wrapped over.
    volatile unsigned int *spin = &mins_idx[gridDim.x];
    while(*spin);

    dist = 999999;  // "infinity", eh.
    hamming = 0;
    for (int i = 0; i < gridDim.x; i++)
    {
	if (mins[i] < dist)
	{
		dist = mins[i];
		hamming = mins_idx[i] + 1; // Note: +1 for matlab
	}
    }

    perm[t] = hamming;

}

