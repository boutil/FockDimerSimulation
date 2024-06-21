package inUtil;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import dimerSim.MarkovSim;
import dimerSim.MarkovSimHex;
import dimerSim.MarkovSimZ2;
import lattices.HexLattice;
import lattices.Visualization;
import lattices.Z2Lattice;



public class HexLatticeTest {
    
    public static void main(String[] args) {
        MarkovSim sim;

        // SchottkyDimersQuad schottkyDimers = buildSchottkyDimers();
        // Z2LatticeFock lattice = new Z2LatticeFock(schottkyDimers, 200, 200);

        HexLattice lattice = new HexLattice(300, 300);


        // // System.out.println(Arrays.deepToString(lattice.faceWeights));

        sim = new MarkovSimHex(lattice);

        sim = loadSim("experimentExport/Hexagon/hexagon300UniformConverged.ser");

        sim.simulate(100000);


        // saveSim(sim, "mySimG2.ser");

        // System.out.println(Arrays.deepToString(lattice.flipFaceWeights));

        Visualization vis = new Visualization(sim);

        // vis.visualizeSim(5);

        // vis.visualizeWeights();

        // vis.visualizeDiscreteAbelMap();

        // vis.visualizeThetaCrossRatio(-3.1691590245331893, 5.071383906313774);

        try {
            System.out.println("simulation done. Saving sim now.");
            vis.saveDimerConfPic("experimentExport/Hexagon/hexagon300UniformConverged.png");
            saveSim(sim, "experimentExport/Hexagon/hexagon300UniformConverged123.ser");
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
    }

    public static void saveSim(MarkovSim sim, String fileName) {
        try {
            ObjectOutputStream out;
            out = new ObjectOutputStream(new FileOutputStream(fileName));
            out.writeObject(sim.lattice);
            out.writeObject(sim.faceStates);
            out.writeObject(sim.insideBoundary);
            out.flush();
            out.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static MarkovSimHex loadSim(String fileName) {
        MarkovSimHex sim = null;
        try {
            ObjectInputStream in;
            in = new ObjectInputStream(new FileInputStream(fileName));
            HexLattice lattice = (HexLattice) in.readObject();
            byte[][] faceStates = (byte[][]) in.readObject();
            boolean[][] insideBoundary = (boolean[][]) in.readObject();
            sim = new MarkovSimHex(lattice, faceStates, insideBoundary);
            in.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return sim;
    }

}
