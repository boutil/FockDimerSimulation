extern "C"
__global__ void flipFaces(int *faceStates, float *faceWeights, int *insideBoundary, float *rand, int size)
{
    int i = (blockIdx.x * blockDim.x + threadIdx.x) / size;
    int j = (blockIdx.x * blockDim.x + threadIdx.x) % size;
    int index = i*size + j;


    if((index < 3*size*size) && (insideBoundary[index] == 1)) {
        float prop = 0;
        if (faceStates[index] == 42) {
            prop = faceWeights[index];
        } else if (faceStates[index] == 21) {
            prop = 1/faceWeights[index];
        }
        if ( rand[index] < prop * 0.9999 ) {
            if ( faceStates[index] == 21 ) { 
                faceStates[index] = 42;
            }
            else if ( faceStates[index] == 42 ) { 
                faceStates[index] = 21;
            }
        }
    }
}

// otherFaceStates 1 and 2 are those of parity one and two larger than theseFaceStates.
// We assume a 3Nx3N grid.
// We use (f(parity) % 3) / 2 to produce a number that's 1 for only one of the parities.
extern "C"
__global__ void consolidateFaces(int *theseFaceStates, int *otherFaceStates1, int *otherFaceStates2, int size, int parity)
{
    int i = ((blockIdx.x * blockDim.x + threadIdx.x) / size) + 1;
    int j = ((blockIdx.x * blockDim.x + threadIdx.x) % size) + 1;
    int index = i*size + j;
    
    if(i < 3 * size - 1 && j < size - 1) {
        theseFaceStates[index] = (otherFaceStates1[(i-1)*size+j] & 1) * 8                                    // W
                                    + (otherFaceStates2[(i+1)*size+j] & 8) / 8                               // E
                                    + (otherFaceStates2[i*size+j - ((parity + i + 2) % 3) / 2] & 2) * 8      // SW
                                    + (otherFaceStates1[i*size+j + ((parity + i) % 3) / 2] & 16) / 8         // NE
                                    + (otherFaceStates1[(i+1)*size+j - ((parity + i + 2) % 3) / 2] & 4) * 8  // SE
                                    + (otherFaceStates2[(i-1)*size+j + ((parity + i) % 3) / 2] & 32) / 8;    // NW
    }
}