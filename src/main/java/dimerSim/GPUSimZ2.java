package dimerSim;

import static jcuda.driver.JCudaDriver.cuLaunchKernel;
import jcuda.Pointer;

import lattices.Lattice;

public class GPUSimZ2 extends GPUSim{

    // We assume a 2Nx2N grid.

    public GPUSimZ2(Lattice lattice, byte[][] faceStates, boolean[][] insideBoundary) {
        super(lattice, faceStates, insideBoundary);
        maxParity = 2;
        ptxFileName = "src/main/java/inUtil/JCudaQuadWalkKernel.ptx";
    }

    @Override
    public void defineKernelParams(int parity) {
        kernelParametersFlip[parity] = Pointer.to(
            Pointer.to(faceStatePointers[parity]),
            Pointer.to(faceWeightPointers[parity]),
            Pointer.to(insideBoundaryPointers[parity]),
            Pointer.to(deviceRandArr),
            Pointer.to(new int[]{N/2})
        );
        kernelParametersCons[parity] = Pointer.to(
            Pointer.to(faceStatePointers[parity]),
            Pointer.to(faceStatePointers[(parity+1) % maxParity]),
            Pointer.to(new int[]{N/2}),
            Pointer.to(new int[]{parity})
        );
    }

    @Override
    public void launchKernels(int parity, int blockSizeX, int gridSizeX) {
        // Call the kernel function.
        cuLaunchKernel(flipFacesF,
            gridSizeX,  1, 1,         // Grid dimension
            blockSizeX, 1, 1,         // Block dimension
            0, null,   // Shared memory size and stream
            kernelParametersFlip[parity], null // Kernel- and extra parameters
        );
        
        cuLaunchKernel(consolidateFacesF,
            gridSizeX,  1, 1,         // Grid dimension
            blockSizeX, 1, 1,         // Block dimension
            0, null,   // Shared memory size and stream
            kernelParametersCons[(parity + 1) % maxParity], null // Kernel- and extra parameters
        );
    }
}
