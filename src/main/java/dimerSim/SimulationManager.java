package dimerSim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import de.jtem.riemann.schottky.SchottkyData;
import de.jtem.riemann.schottky.SchottkyDimers;
import de.jtem.riemann.schottky.SchottkyDimersDoubleCoverUnitary;
import de.jtem.riemann.schottky.SchottkyDimersQuad;
import de.jtem.riemann.schottky.SchottkyDimersUnitary;
import lattices.HexLattice;
import lattices.Lattice;
import lattices.Visualization;
import lattices.Z2Lattice;

public class SimulationManager {
    // class responsible for launching and managing simulations.
    // Runs simulations corresponding to a single schottkyDimers.
    // Folder structure:
    // | schottky.ser
    // | 300x300
    // --| evolutionPics
    // ----| PicturesOverTime
    // --| runDate1
    // --| runDate2
    // --| .....
    // | 500x500
    // --| runDate1
    // --| runDate2
    // --| .....

    private SchottkyDimers schottkyDimers;
    private MarkovSim sim;
    private Visualization vis;

    private String baseFolder;

    // How many simulation steps between saving a picture.
    private boolean saveEvolutionPics;
    private int savePictureInterval;


    Map<Integer, String> hexagonUnifPaths = new HashMap<Integer, String>();


    private final String Hexagon300Unif = "";
    private final String Hexagon500Unif = "";


    public SimulationManager(String folder) {
        baseFolder = folder;
        schottkyDimers = loadSchottky();
        init();
    }

    public SimulationManager(SchottkyDimers schottky, String folder) {
        baseFolder = folder;
        new File(baseFolder).mkdirs();
        schottkyDimers = schottky;
        saveSchottky(schottky, new File(folder, "schottky.ser"));
        init();
    }

    private void init() {
        vis = new Visualization(sim, schottkyDimers);
        hexagonUnifPaths.put(300, Hexagon300Unif);
        hexagonUnifPaths.put(500, Hexagon500Unif);

    }

    public void setSavePictureInterval(int interval) {
        // stop saving interval by interval = 0
        savePictureInterval = interval;
        saveEvolutionPics = interval != 0;
    }

    private MarkovSim createLattice(int size) {
        // Need a more versatile way of constructing lattices here.
        return loadSim(hexagonUnifPaths.get(size));
    }

    public void simulateAndSave(int numSteps, int size) {
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy-MM-dd-HH-mm-ss");
        LocalDateTime now = LocalDateTime.now();
        File sizeFolder = new File(baseFolder, size + "x" + size);
        File simFolder = new File(sizeFolder, dtf.format(now));
        File evoFolder = new File(sizeFolder, "evolutionPics");
        // create new folder for 
        simFolder.mkdirs();
        evoFolder.mkdirs();

        sim = loadNewestSim(sizeFolder, size);
        if(!saveEvolutionPics) {
            sim.simulate(numSteps);
        } else {
            int numFilesInEvo = evoFolder.listFiles().length;
            for (int i = 0; i < numSteps/savePictureInterval; i++) {
                sim.simulate(savePictureInterval);
                String filePath = new File(evoFolder, (i + numFilesInEvo) + ".png").getPath();
                try {
                    vis.saveDimerConfPic(filePath, true);
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
        }
        saveSim(simFolder);
    }
    
    private void saveSim(File simFolder) {
        try {
            String info = "[" + sim.lattice.N + "x" + sim.lattice.M + "]";
            saveSim(sim, simFolder + "/sim.ser");
            vis.saveDimerConfPic(simFolder + "/dimerConf.png", false);
            vis.saveWeightsPic(simFolder + "/weights.png");
            vis.saveDimerConfPic(simFolder + "/dimerConfPred.png", true);
            // vis.saveAmoebaPic(simFolder + "/amoebaPic" + info + ".png");
            // vis.saveAztecPic(schottkyDimers, simFolder + "/aztecPic" + info + ".png");
        } catch (IOException e) {
            // TODO: handle exception
        }
    }
    
    
    private MarkovSim loadNewestSim(File sizeDir, int size) {
        if(!sizeDir.mkdirs()) {
            // directory exists so load newest sim.
            List<File> simDirs = Arrays.asList(sizeDir.listFiles());
            Collections.sort(simDirs);
            File newest = simDirs.get(simDirs.size());
            return loadSim(newest.getPath() + "sim.ser");   
        }
        return createLattice(size);
    }

    

    private SchottkyDimers loadSchottky() {
        SchottkyDimers schottkyDimers = null;
        try {
            ObjectInputStream in;
            in = new ObjectInputStream(new FileInputStream(baseFolder + "schottky.ser"));
            double[] uniformizationData = (double[]) in.readObject();
            double[][] angles = (double[][]) in.readObject();
            if(angles.length == 3) {
                schottkyDimers = new SchottkyDimersUnitary(new SchottkyData(uniformizationData), angles);
            } else if (angles.length == 4) {
                schottkyDimers = new SchottkyDimersQuad(new SchottkyData(uniformizationData), angles);
            } else if (angles.length == 6) {
                schottkyDimers = new SchottkyDimersDoubleCoverUnitary(new SchottkyData(uniformizationData), angles);
            }
            in.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return schottkyDimers;
    }

    private void saveSchottky(SchottkyDimers schottkyDimers, File file) {
        try {
            ObjectOutputStream out;
            out = new ObjectOutputStream(new FileOutputStream(file));
            out.writeObject(schottkyDimers.getUniformizationData());
            out.writeObject(schottkyDimers.angles);
            out.flush();
            out.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private void saveSim(MarkovSim sim, String fileName) {
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

    private MarkovSim loadSim(String fileName) {
        MarkovSim sim = null;
        try {
            ObjectInputStream in;
            in = new ObjectInputStream(new FileInputStream(fileName));
            Lattice l = (Lattice) in.readObject();
            byte[][] faceStates = (byte[][]) in.readObject();
            boolean[][] insideBoundary = (boolean[][]) in.readObject();
            if (l.getClass().isAssignableFrom(HexLattice.class)) {
                HexLattice lattice = (HexLattice) l;
                sim = new MarkovSimHex(lattice, faceStates, insideBoundary);
            } else if (l.getClass().isAssignableFrom(Z2Lattice.class)) {
                Z2Lattice lattice = (Z2Lattice) l;
                sim = new MarkovSimZ2(lattice, faceStates, insideBoundary);
            }
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
