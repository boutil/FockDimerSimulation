package dimerSim;

import java.util.Random;

import static jcuda.jcurand.JCurand.curandCreateGenerator;
import static jcuda.jcurand.JCurand.curandDestroyGenerator;
import static jcuda.jcurand.JCurand.curandGenerateUniform;
import static jcuda.jcurand.JCurand.curandSetPseudoRandomGeneratorSeed;
import static jcuda.jcurand.curandRngType.CURAND_RNG_PSEUDO_DEFAULT;
import static jcuda.runtime.JCuda.cudaFree;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpy;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;

import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUdeviceptr;
import jcuda.jcurand.JCurand;
import jcuda.runtime.JCuda;
import jcuda.jcurand.curandGenerator;


import lattices.Lattice;

public class GPUSim {
    // class that will run simulation steps on CUDA compatible GPU.

    Lattice lattice;

    byte[][] faceStates;
    boolean[][] insideBoundary;
    
    Random rand;


    // Cuda
    int evenSize;
    int oddSize;

    CUdeviceptr faceWeightsEvenPointer;
    CUdeviceptr faceStatesEvenPointer;
    CUdeviceptr insideBoundaryEvenPointer;
    CUdeviceptr faceWeightsOddPointer;
    CUdeviceptr faceStatesOddPointer;
    CUdeviceptr insideBoundaryOddPointer;

    Pointer deviceRandArr;
    curandGenerator generator;


    public GPUSim(Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        this.lattice = lattice;
        this.faceStates = faceStates;
        this.insideBoundary = insideBoundary;
        rand = new Random();
    }

    private void initCuda() {
        // Copy memory and load kernel functions
        int N = lattice.N;

        // separate 2d arrays into odd and even arrays.
        // boolean does not exist on GPU -> transform to byte. -> could also encode that into the faceState. Maybe better?
        evenSize = (N * N) / 2;
        oddSize = (N * N) / 2 + (N % 2);
        float[] faceWeightsEven = new float[evenSize];
        byte[] insideBoundaryEven = new byte[evenSize];
        byte[] faceStatesEven = new byte[evenSize];
        float[] faceWeightsOdd = new float[oddSize];
        byte[] insideBoundaryOdd = new byte[oddSize];
        byte[] faceStatesOdd = new byte[oddSize];
        int currentIndex = 0;
        for (int i = 0; i < lattice.flipFaceWeights.length; i++) {
            for (int j = 0; j < lattice.flipFaceWeights[i].length; j++) {
                if((i + j) % 2 == 0) {
                    faceWeightsEven[currentIndex] = (float)lattice.flipFaceWeights[i][j];
                    faceStatesEven[currentIndex] = faceStates[i][j];
                    if(insideBoundary[i][j]) {
                        insideBoundaryEven[currentIndex] = 0b1;
                    } else {
                        insideBoundaryEven[currentIndex] = 0b0;
                    }
                } else {
                    faceWeightsOdd[currentIndex] = (float)lattice.flipFaceWeights[i][j];
                    faceStatesOdd[currentIndex] = faceStates[i][j];
                    if(insideBoundary[i][j]) {
                        insideBoundaryOdd[currentIndex++] = 0b1;
                    } else {
                        insideBoundaryOdd[currentIndex++] = 0b0;
                    }
                }
            }
        }

        // Create on device pointers
        faceWeightsEvenPointer = new CUdeviceptr();
        faceStatesEvenPointer = new CUdeviceptr();
        insideBoundaryEvenPointer = new CUdeviceptr();
        faceWeightsOddPointer = new CUdeviceptr();
        faceStatesOddPointer = new CUdeviceptr();
        insideBoundaryOddPointer = new CUdeviceptr();

        // Copy data to GPU
        cudaMalloc(faceWeightsEvenPointer, evenSize * Sizeof.FLOAT);
        cudaMemcpy(faceWeightsEvenPointer, Pointer.to(faceWeightsEven), evenSize * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
        cudaMalloc(faceStatesEvenPointer, evenSize * Sizeof.BYTE);
        cudaMemcpy(faceStatesEvenPointer, Pointer.to(faceStatesEven), evenSize * Sizeof.BYTE, cudaMemcpyDeviceToHost);
        cudaMalloc(insideBoundaryEvenPointer, evenSize * Sizeof.BYTE);
        cudaMemcpy(insideBoundaryEvenPointer, Pointer.to(insideBoundaryEven), evenSize * Sizeof.BYTE, cudaMemcpyDeviceToHost);

        cudaMalloc(faceWeightsOddPointer, oddSize * Sizeof.FLOAT);
        cudaMemcpy(faceWeightsOddPointer, Pointer.to(faceWeightsOdd), oddSize * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
        cudaMalloc(faceStatesOddPointer, oddSize * Sizeof.BYTE);
        cudaMemcpy(faceStatesOddPointer, Pointer.to(faceStatesOdd), oddSize * Sizeof.BYTE, cudaMemcpyDeviceToHost);
        cudaMalloc(insideBoundaryOddPointer, oddSize * Sizeof.BYTE);
        cudaMemcpy(insideBoundaryOddPointer, Pointer.to(insideBoundaryOdd), oddSize * Sizeof.BYTE, cudaMemcpyDeviceToHost);


        
        // init random generator;
        generator = new curandGenerator();

        deviceRandArr = new Pointer();
        cudaMalloc(deviceRandArr, oddSize * Sizeof.FLOAT);
        curandCreateGenerator(generator, CURAND_RNG_PSEUDO_DEFAULT);
        curandSetPseudoRandomGeneratorSeed(generator, 1234);
        // To generate random numbers:
        // curandGenerateUniform(generator, deviceRandArr, N);
    }

    private void cleanUpCuda() {
        curandDestroyGenerator(generator);
        cudaFree(deviceRandArr);
    }



    public void simulate(int steps) {
        initCuda();

        for (int i = 0; i < steps; i++) {
            int parity = rand.nextInt(2);
            if (parity == 0) {
                curandGenerateUniform(generator, deviceRandArr, evenSize);
                // cuda.flipFaces(even)
                // cuda.consolidateFaces(odd)
            } else {
                curandGenerateUniform(generator, deviceRandArr, oddSize);
            }
        }

        cleanUpCuda();
    }

}
