
extern "C"
__global__ void flipFaces(int *faceStates, float *faceWeights, int *insideBoundary, float *rand, int size)
{
    int i = (blockIdx.x * blockDim.x + threadIdx.x) / size;
    int j = (blockIdx.x * blockDim.x + threadIdx.x) % size;
    int index = i*size + j;
    // faceStates[index] = index;


    if((index < 2*size*size) && (insideBoundary[index] == 1)) {
        float prop = 0;
        if (faceStates[index] == 10) {
            prop = faceWeights[index];
        } else if (faceStates[index] == 5) {
            prop = 1/faceWeights[index];
        }
        if ( rand[index] < prop * 0.9999 ) {
            if ( faceStates[index] == 5 ) { 
                faceStates[index] = 10;
            }
            else if ( faceStates[index] == 10 ) { 
                faceStates[index] = 5;
            }
        }
    }
}


extern "C"
__global__ void consolidateFaces(int *theseFaceStates, int *otherFaceStates, int size, int parity)
{
    int i = ((blockIdx.x * blockDim.x + threadIdx.x) / size) + 1;
    int j = ((blockIdx.x * blockDim.x + threadIdx.x) % size) + 1;
    int index = i*size + j;
    
    if(i < 2 * size - 1 && j < size - 1) {
        theseFaceStates[index] = (otherFaceStates[(i-1)*size+j] & 4)/4
                                    + 4*(otherFaceStates[(i+1)*size+j] & 1) 
                                    + (otherFaceStates[i*size+j - (parity + i + 1) % 2] & 8)/4 
                                    + 4*(otherFaceStates[i*size+j + (parity + i) % 2] & 2);
    }
}