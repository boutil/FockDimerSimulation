package inUtil;

import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

import javax.imageio.ImageIO;

import com.squareup.gifencoder.FloydSteinbergDitherer;
import com.squareup.gifencoder.GifEncoder;
import com.squareup.gifencoder.ImageOptions;

public class GifConvert {
    public static void main(String[] args) throws Exception {
        GifConvert gifTutorial = new GifConvert();
        gifTutorial.createAnimatedGif();
    }

    // faster way to convert to gif:
    // ffmpeg -i %06d.png -filter_complex "fps=5,scale=960:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=32[p];[s1][p]paletteuse=dither=bayer" output.gif
 

    int width;
    int height;
 
    private void createAnimatedGif() throws Exception {
        String folderPath = "experimentExport/Hexagon/03_convergenceToUniform/300x300/evolutionPics";
        File folder = new File(folderPath);
        List<File> simDirs = Arrays.asList(folder.listFiles());
        Collections.sort(simDirs);

        //The GIF image will be created with file name "my_animated_image.gif"
        try (FileOutputStream outputStream = new FileOutputStream(new File(folderPath, "animation.gif"))) {
            ImageOptions options = new ImageOptions();

            //Set delay between each frame
            options.setDelay(100, TimeUnit.MILLISECONDS);
            //Use Floyd Steinberg dithering as it yields the best quality
            options.setDitherer(FloydSteinbergDitherer.INSTANCE);

            // BufferedImage bI = ImageIO.read(simDirs.get(0));
            // width = bI.getWidth();
            // height = bI.getHeight();
            width = 540 * 2;
            height = 311 * 2;

            //Create GIF encoder with same dimension as of the source images
            GifEncoder enc = new GifEncoder(outputStream, width, height, 0);

            int i = 0;
            for (File file : simDirs) {
                System.out.println(i++);
                enc.addImage(convertImageToArray(file), options);
            }
            enc.finishEncoding(); //Start the encoding
        }
    }
 
  /**
   * Convert BufferedImage into RGB pixel array
   */
    public int[][] convertImageToArray(File file) throws IOException {
        BufferedImage bufferedImage = ImageIO.read(file);
        // resizing:
        BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = (Graphics2D)bi.createGraphics();
        g2d.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING,
                RenderingHints.VALUE_RENDER_QUALITY));
        boolean b = g2d.drawImage(bufferedImage, 0, 0, width, height, null);

        int[][] rgbArray = new int[bi.getHeight()][bi.getWidth()];
        for (int i = 0; i < bi.getHeight(); i++) {
            for (int j = 0; j < bi.getWidth(); j++) {
                
                rgbArray[i][j] = bi.getRGB(j, i);
            }
        }
        return rgbArray;
    }
}
