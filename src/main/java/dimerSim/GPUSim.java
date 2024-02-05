package dimerSim;

import java.util.Random;

import org.jzy3d.plot3d.pipelines.NotImplementedException;

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
import lattices.Lattice;
import jcuda.jcurand.curandGenerator;

public class GPUSim {
    
    private Lattice lattice;

    private byte[][] faceStates;
    private boolean[][] insideBoundary;
    
    private Random rand;


    // Cuda
    int N;
    int arrSize;
    int maxParity;

    public String ptxFileName;

    protected CUdeviceptr[] faceWeightPointers;
    protected CUdeviceptr[] faceStatePointers;
    protected CUdeviceptr[] insideBoundaryPointers;
 
    private CUcontext context;
    private CUmodule walkModule;
 
    protected CUfunction flipFacesF;
    protected CUfunction consolidateFacesF;
 
    protected Pointer[] kernelParametersFlip;
    protected Pointer[] kernelParametersCons;
     
    protected Pointer deviceRandArr;
    private curandGenerator generator;

    // We assume a maxParity x maxParity grid.

    public GPUSim(Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        this.lattice = lattice;
        this.faceStates = faceStates;
        this.insideBoundary = insideBoundary;
        rand = new Random();
        N = lattice.N;
    }


    // For compiling shader in Windows: In developer terminal:
    // "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64
    // nvcc -ptx JCudaQuadWalkKernel.cu -o JCudaQuadWalkKernel.ptx


    public void simulate(int steps) {
        try {

            initCuda();
            
            float[] floatArrToShow = new float[arrSize];
            
            int blockSizeX = 256;
            int gridSizeX = (int)Math.ceil((double)arrSize / blockSizeX);
            
            
            int reportFreq = steps/10;
            long time = System.currentTimeMillis();
            for (int i = 0; i < steps; i++) {
                if ((i % reportFreq == 0)) {
                    double timeForAvg1000Steps = (double)(System.currentTimeMillis() - time) * 1000 / reportFreq;
                    System.out.println("Done with " + i + " steps." + " Average time per 1000 markovSteps: " + Math.ceil(timeForAvg1000Steps) + ". Time left: " + Math.ceil(((steps - i) * timeForAvg1000Steps / 1000000)) + " seconds.");
                    time = System.currentTimeMillis();
                }
                
                int parity = rand.nextInt(maxParity);
                curandSetGeneratorOffset(generator, ((long) i * (long) arrSize) % Long.MAX_VALUE);
                
                curandGenerateUniform(generator, deviceRandArr, arrSize);
                
                launchKernels(parity, blockSizeX, gridSizeX);
                
            }
            // Copy data back to faceStates on host.
            copyDataBack();
        } finally {
            cleanUpCuda();
        }
    }

    protected void defineKernelParams(int parity) {
        throw new NotImplementedException();
    }

    protected void launchKernels(int parity, int blocksizeX, int gridSizeX) {
        throw new NotImplementedException();
    }



    private void copyDataBack() {
        int[][] faceStatesArr = new int[maxParity][arrSize];
        for (int i = 0; i < maxParity; i++) {
            cudaMemcpy(Pointer.to(faceStatesArr[i]), faceStatePointers[i], arrSize * Sizeof.INT, cudaMemcpyDeviceToHost);
        }
        for (int i = 0; i < faceStates.length; i++) {
            for (int j = 0; j < faceStates[i].length; j++) {
                int parity = Math.floorMod(j - i, maxParity);
                faceStates[i][j] = (byte)faceStatesArr[parity][(i * N + j) / maxParity];
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

        // separate 2d arrays into odd and even 1d arrays.
        // boolean does not exist on GPU -> transform to int. -> could also encode that into one bit of the faceState. Maybe better?
        arrSize = (N * N) / maxParity;
        float[][] faceWeightsArr = new float[maxParity][];
        int[][] faceStatesArr = new int[maxParity][];
        int[][] insideBdryArr = new int[maxParity][];
        for (int i = 0; i < maxParity; i++) {
            faceWeightsArr[i] = new float[arrSize];
            faceStatesArr[i] = new int[arrSize];
            insideBdryArr[i] = new int[arrSize];
        }
        // May want to pad with zeros to achieve blocksize*gridsize number of elements exactly. Then kernel doesn't need to check for i < size
        for (int i = 0; i < lattice.flipFaceWeights.length; i++) {
            for (int j = 0; j < lattice.flipFaceWeights[i].length; j++) {
                int parity = Math.floorMod(j - i, maxParity);
                faceWeightsArr[parity][(i * N + j) / maxParity] = (float)lattice.flipFaceWeights[i][j];
                faceStatesArr[parity][(i * N + j) / maxParity] = faceStates[i][j];
                if(insideBoundary[i][j]) {
                    insideBdryArr[parity][(i * N + j) / maxParity] = 1;
                } else {
                    insideBdryArr[parity][(i * N + j) / maxParity] = 0;
                }
            }
        }

        faceWeightPointers = new CUdeviceptr[maxParity];
        faceStatePointers = new CUdeviceptr[maxParity];
        insideBoundaryPointers = new CUdeviceptr[maxParity];
        // Create on device pointers
        for (int i = 0; i < maxParity; i++) {
            faceWeightPointers[i] = new CUdeviceptr();
            faceStatePointers[i] = new CUdeviceptr();
            insideBoundaryPointers[i] = new CUdeviceptr();
        }

        // Copy data to GPU
        for (int i = 0; i < maxParity; i++) {
            cudaMalloc(faceWeightPointers[i], arrSize * Sizeof.FLOAT);
            cudaMemcpy(faceWeightPointers[i], Pointer.to(faceWeightsArr[i]), arrSize * Sizeof.FLOAT, cudaMemcpyHostToDevice);
            cudaMalloc(faceStatePointers[i], arrSize * Sizeof.INT);
            cudaMemcpy(faceStatePointers[i], Pointer.to(faceStatesArr[i]), arrSize * Sizeof.INT, cudaMemcpyHostToDevice);
            cudaMalloc(insideBoundaryPointers[i], arrSize * Sizeof.INT);
            cudaMemcpy(insideBoundaryPointers[i], Pointer.to(insideBdryArr[i]), arrSize * Sizeof.INT, cudaMemcpyHostToDevice);
        }
        

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
        cudaMalloc(deviceRandArr, arrSize * Sizeof.FLOAT);
        curandCreateGenerator(generator, CURAND_RNG_PSEUDO_DEFAULT);
        curandSetPseudoRandomGeneratorSeed(generator, 1234);

        kernelParametersFlip = new Pointer[maxParity];
        kernelParametersCons = new Pointer[maxParity];
        
        for (int i = 0; i < maxParity; i++) {
            defineKernelParams(i);
        }
        // To generate random numbers:
        // curandGenerateUniform(generator, deviceRandArr, N);
    }

    private void cleanUpCuda() {
        cuModuleUnload(walkModule);

        for (int i = 0; i < maxParity; i++) {
            cuMemFree(faceWeightPointers[i]);
            cuMemFree(faceStatePointers[i]);
            cuMemFree(insideBoundaryPointers[i]);
        }

        curandDestroyGenerator(generator);
        cudaFree(deviceRandArr);

        cuCtxDestroy(context);
    }
}
