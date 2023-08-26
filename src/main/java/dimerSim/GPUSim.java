package dimerSim;

import java.util.Random;

import static jcuda.jcurand.JCurand.curandCreateGenerator;
import static jcuda.jcurand.JCurand.curandDestroyGenerator;
import static jcuda.jcurand.JCurand.curandGenerateUniform;
import static jcuda.jcurand.JCurand.curandSetGeneratorOffset;
import static jcuda.jcurand.JCurand.curandSetPseudoRandomGeneratorSeed;
import static jcuda.jcurand.curandRngType.CURAND_RNG_PSEUDO_DEFAULT;
import static jcuda.runtime.JCuda.cudaFree;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpy;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;
import static jcuda.driver.JCudaDriver.cuMemFree;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleLoad;
import static jcuda.driver.JCudaDriver.cuModuleUnload;
import static jcuda.driver.JCudaDriver.cuCtxSynchronize;
import static jcuda.driver.JCudaDriver.cuLaunchKernel;
import static jcuda.driver.JCudaDriver.cuDeviceGet;
import static jcuda.driver.JCudaDriver.cuInit;
import static jcuda.driver.JCudaDriver.cuCtxCreate;
import static jcuda.driver.JCudaDriver.cuCtxDestroy;
import jcuda.driver.CUmodule;
import jcuda.driver.CUfunction;




import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.jcurand.JCurand;
import jcuda.runtime.JCuda;
import jcuda.jcurand.curandGenerator;


import lattices.Lattice;

public class GPUSim {
    // class that will run simulation steps on CUDA compatible GPU.



    // For compiling shader in Windows: In developer terminal:
    // "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64 
    // nvcc -ptx JCudaReductionKernel.cu -o JCudaReductionKernel.ptx



    private Lattice lattice;

    private byte[][] faceStates;
    private boolean[][] insideBoundary;
    
    private Random rand;


    // Cuda
    int N;
    int evenSize;
    int oddSize;

    private CUdeviceptr faceWeightsEvenPointer;
    private CUdeviceptr faceStatesEvenPointer;
    private CUdeviceptr insideBoundaryEvenPointer;
    private CUdeviceptr faceWeightsOddPointer;
    private CUdeviceptr faceStatesOddPointer;
    private CUdeviceptr insideBoundaryOddPointer;
 
    private CUcontext context;
    private CUmodule walkModule;
 
    private CUfunction flipFacesF;
    private CUfunction consolidateFacesF;
 
    private Pointer kernelParametersFlipEven;
    private Pointer kernelParametersFlipOdd;
     
    private Pointer kernelParametersConsEven;
    private Pointer kernelParametersConsOdd;
 
    private Pointer deviceRandArr;
    private curandGenerator generator;
    
    
    public GPUSim(Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        this.lattice = lattice;
        this.faceStates = faceStates;
        this.insideBoundary = insideBoundary;
        rand = new Random();
    }


    public void simulate(int steps) {
        initCuda();

        int[] intArrToShow = new int[evenSize];
        float[] floatArrToShow = new float[evenSize];
        
        cudaMemcpy(Pointer.to(intArrToShow), faceStatesEvenPointer, evenSize * Sizeof.INT, cudaMemcpyDeviceToHost);
        intArrToShow = new int[evenSize];
        
        int blockSizeX = 256;
        int gridSizeX = (int)Math.ceil((double)evenSize / blockSizeX);
        

        int reportFreq = steps/10;
        long time = System.currentTimeMillis();
        for (int i = 0; i < steps; i++) {
            if ((i % reportFreq == 0)) {
                double timeForAvg1000Steps = Math.max(((double)(System.currentTimeMillis() - time)) * 1000 / reportFreq, 1000);
                System.out.println("Done with " + i + " steps." + " Average time per 1000 markovSteps: " + (int) timeForAvg1000Steps + ". Time left: " + (int)((steps - i) * timeForAvg1000Steps / 1000000) + " seconds.");
                time = System.currentTimeMillis();
            }

            int parity = rand.nextInt(2);
            curandSetGeneratorOffset(generator, ((long) i * (long) evenSize) % Long.MAX_VALUE);
            if (parity == 0) {
                curandGenerateUniform(generator, deviceRandArr, evenSize);
                
                // Call the kernel function.
                cuLaunchKernel(flipFacesF,
                gridSizeX,  1, 1,         // Grid dimension
                blockSizeX, 1, 1,         // Block dimension
                0, null,   // Shared memory size and stream
                kernelParametersFlipEven, null // Kernel- and extra parameters
                );
                // cuCtxSynchronize();
                
                cuLaunchKernel(consolidateFacesF,
                    gridSizeX,  1, 1,         // Grid dimension
                    blockSizeX, 1, 1,         // Block dimension
                    0, null,   // Shared memory size and stream
                    kernelParametersConsOdd, null // Kernel- and extra parameters
                );
                // cuCtxSynchronize();

            } else {
                curandGenerateUniform(generator, deviceRandArr, oddSize);
                // Call the kernel function.
                cuLaunchKernel(flipFacesF,
                    gridSizeX,  1, 1,         // Grid dimension
                    blockSizeX, 1, 1,         // Block dimension
                    0, null,   // Shared memory size and stream
                    kernelParametersFlipOdd, null // Kernel- and extra parameters
                );

                // cuCtxSynchronize();
                cuLaunchKernel(consolidateFacesF,
                    gridSizeX,  1, 1,         // Grid dimension
                    blockSizeX, 1, 1,         // Block dimension
                    0, null,   // Shared memory size and stream
                    kernelParametersConsEven, null // Kernel- and extra parameters
                );
                // cuCtxSynchronize();
            }
        }

        // Copy data back to faceStates on host.
        copyDataBack();
        

        cleanUpCuda();
    }

    public byte[][] getFaceStates() {
        return faceStates;
    }


    private void copyDataBack() {
        int[] faceStatesEven = new int[evenSize];
        int[] faceStatesOdd = new int[oddSize];
        cudaMemcpy(Pointer.to(faceStatesEven), faceStatesEvenPointer, evenSize * Sizeof.INT, cudaMemcpyDeviceToHost);
        cudaMemcpy(Pointer.to(faceStatesOdd), faceStatesOddPointer, oddSize * Sizeof.INT, cudaMemcpyDeviceToHost);
        for (int i = 0; i < faceStates.length; i++) {
            for (int j = 0; j < faceStates[i].length; j++) {
                if((i+j) % 2 ==0) {
                    faceStates[i][j] = (byte)faceStatesEven[(i * N + j) / 2];
                } else {
                    faceStates[i][j] = (byte)faceStatesOdd[(i * N + j) / 2];
                }
            }
        }
    }


    private void initCuda() {
        // Copy memory and load kernel functions
        cuInit(0);
        CUdevice device = new CUdevice();
        cuDeviceGet(device, 0);
        context = new CUcontext();
        cuCtxCreate(context, 0, device);

        N = lattice.N;

        // separate 2d arrays into odd and even 1d arrays.
        // boolean does not exist on GPU -> transform to int. -> could also encode that into one bit of the faceState. Maybe better?
        evenSize = (N * N) / 2 + (N % 2);
        oddSize = (N * N) / 2;
        float[] faceWeightsEven = new float[evenSize];
        int[] insideBoundaryEven = new int[evenSize];
        int[] faceStatesEven = new int[evenSize];
        float[] faceWeightsOdd = new float[oddSize];
        int[] insideBoundaryOdd = new int[oddSize];
        int[] faceStatesOdd = new int[oddSize];
        // int currentIndex = 0;
        for (int i = 0; i < lattice.flipFaceWeights.length; i++) {
            for (int j = 0; j < lattice.flipFaceWeights[i].length; j++) {
                if((i + j) % 2 == 0) {
                    faceWeightsEven[(i * N + j) / 2] = (float)lattice.flipFaceWeights[i][j];
                    faceStatesEven[(i * N + j) / 2] = faceStates[i][j];
                    // insideBoundaryEven[(i * N + j) / 2] = insideBoundary[i][j];
                    if(insideBoundary[i][j]) {
                        insideBoundaryEven[(i * N + j) / 2] = 1;
                    } else {
                        insideBoundaryEven[(i * N + j) / 2] = 0;
                    }
                } else {
                    faceWeightsOdd[(i * N + j) / 2] = (float)lattice.flipFaceWeights[i][j];
                    faceStatesOdd[(i * N + j) / 2] = faceStates[i][j];
                    // insideBoundaryOdd[(i * N + j) / 2] = insideBoundary[i][j];
                    if(insideBoundary[i][j]) {
                        insideBoundaryOdd[(i * N + j) / 2] = 1;
                    } else {
                        insideBoundaryOdd[(i * N + j) / 2] = 0;
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
        cudaMemcpy(faceWeightsEvenPointer, Pointer.to(faceWeightsEven), evenSize * Sizeof.FLOAT, cudaMemcpyHostToDevice);
        cudaMalloc(faceStatesEvenPointer, evenSize * Sizeof.INT);
        cudaMemcpy(faceStatesEvenPointer, Pointer.to(faceStatesEven), evenSize * Sizeof.INT, cudaMemcpyHostToDevice);
        cudaMalloc(insideBoundaryEvenPointer, evenSize * Sizeof.INT);
        cudaMemcpy(insideBoundaryEvenPointer, Pointer.to(insideBoundaryEven), evenSize * Sizeof.INT, cudaMemcpyHostToDevice);
        
        cudaMalloc(faceWeightsOddPointer, oddSize * Sizeof.FLOAT);
        cudaMemcpy(faceWeightsOddPointer, Pointer.to(faceWeightsOdd), oddSize * Sizeof.FLOAT, cudaMemcpyHostToDevice);
        cudaMalloc(faceStatesOddPointer, oddSize * Sizeof.INT);
        cudaMemcpy(faceStatesOddPointer, Pointer.to(faceStatesOdd), oddSize * Sizeof.INT, cudaMemcpyHostToDevice);
        cudaMalloc(insideBoundaryOddPointer, oddSize * Sizeof.INT);
        cudaMemcpy(insideBoundaryOddPointer, Pointer.to(insideBoundaryOdd), oddSize * Sizeof.INT, cudaMemcpyHostToDevice);
        
        String ptxFileName = "src/main/java/inUtil/JCudaQuadWalkKernel.ptx";
        walkModule = new CUmodule();
        cuModuleLoad(walkModule, ptxFileName);
        
        // Obtain a function pointer to the "reduce" function.
        flipFacesF = new CUfunction();
        cuModuleGetFunction(flipFacesF, walkModule, "flipFaces");
        
        consolidateFacesF = new CUfunction();
        cuModuleGetFunction(consolidateFacesF, walkModule, "consolidateFaces");


        // init random generator;
        generator = new curandGenerator();

        deviceRandArr = new CUdeviceptr();
        cudaMalloc(deviceRandArr, oddSize * Sizeof.FLOAT);
        curandCreateGenerator(generator, CURAND_RNG_PSEUDO_DEFAULT);
        curandSetPseudoRandomGeneratorSeed(generator, 1234);


        kernelParametersFlipEven = Pointer.to(
            Pointer.to(faceStatesEvenPointer),
            Pointer.to(faceWeightsEvenPointer),
            Pointer.to(insideBoundaryEvenPointer),
            Pointer.to(deviceRandArr),
            Pointer.to(new int[]{N/2})
            );
            
        kernelParametersFlipOdd = Pointer.to(
            Pointer.to(faceStatesOddPointer),
            Pointer.to(faceWeightsOddPointer),
            Pointer.to(insideBoundaryOddPointer),
            Pointer.to(deviceRandArr),
            Pointer.to(new int[]{N/2})
        );

        kernelParametersConsEven = Pointer.to(
            Pointer.to(faceStatesEvenPointer),
            Pointer.to(faceStatesOddPointer),
            Pointer.to(new int[]{N/2}),
            Pointer.to(new int[]{0})
        );

        kernelParametersConsOdd = Pointer.to(
            Pointer.to(faceStatesOddPointer),
            Pointer.to(faceStatesEvenPointer),
            Pointer.to(new int[]{N/2}),
            Pointer.to(new int[]{1})
        );
        
        // To generate random numbers:
        // curandGenerateUniform(generator, deviceRandArr, N);
    }

    private void cleanUpCuda() {
        cuModuleUnload(walkModule);

        cuMemFree(faceWeightsEvenPointer);
        cuMemFree(faceStatesEvenPointer);
        cuMemFree(insideBoundaryEvenPointer);
        cuMemFree(faceWeightsOddPointer);
        cuMemFree(faceStatesOddPointer);
        cuMemFree(insideBoundaryOddPointer);

        curandDestroyGenerator(generator);
        cudaFree(deviceRandArr);

        cuCtxDestroy(context);
    }
}
