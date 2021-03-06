package phase_field_model;

import javax.swing.*;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.Timer;
import java.util.TimerTask;

/**
 * PottsViewPanel.java
 * 
 * Paint the lattice configuration on screen
 * 
 * @author Michael Chiang
 *
 */

@SuppressWarnings("serial")
public class PhaseFieldViewPanel extends JPanel 
  implements PhaseFieldDataListener {
	
	private PhaseFieldModel model;
	private int numOfColourBins = 20;
	private Color [][] colours;
	private double minVal = 0.0;
	private double maxVal = 1.0;
	private double incVal = (maxVal-minVal)/numOfColourBins;
	private BufferedImage fg = null;
	private Timer timer = null;
	
	
	/**
	 * Initialise the view panel
	 * @param model the spin model to be displayed on screen
	 */
	public PhaseFieldViewPanel(PhaseFieldModel model){
		this.model = model;
		this.model.addDataListener(this);
		setColours();
	}
	
	/**
	 * Set the spin model to be displayed on screen
	 * @param model spin model
	 */
	public void setModel(PhaseFieldModel model){
		this.model.removeDataListener(this);
		this.model = model;
		this.model.addDataListener(this);
		setColours();
		repaint();
	}
	
	/**
	 * Set the colour for each spin
	 */
	public void setColours(){
		colours = new Color [2][numOfColourBins];
		int r, g, b;
		for (int i = 0; i < 2; i++){
		  for (int j = 0; j < numOfColourBins; j++){
		    double period = (i+1)*Math.PI;
		    r = (int) (Math.sin(period/(double) j + 0)*127.0 + 128.0);
		    g = (int) (Math.sin(period/(double) j + period/3.0)*127.0 + 128.0);
		    b = (int) (Math.sin(period/(double) j + 2.0*period/3.0)*127.0 + 128.0);
		    colours[i][j] = new Color(r, g, b);
		  }
		}
	}
	
	/**
	 * Begin painting the model configuration on screen with fixed 
	 * update rate
	 */
	public void initImage(){
		int width = model.getLx();
		int height = model.getLy();
		
		fg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			
		drawImage();
		
		final JPanel panel = this;
		panel.getGraphics().drawImage(fg, 0, 
				panel.getInsets().top, panel.getWidth(), 
				panel.getHeight() - panel.getInsets().top, null);
		
		timer = new Timer();
		timer.scheduleAtFixedRate(new TimerTask() {
			public void run() {
				panel.getGraphics().drawImage(
						fg, 0, panel.getInsets().top, 
						panel.getWidth(), 
						panel.getHeight() - panel.getInsets().top, null);
			}
		}, 0, 33);
	}
	
	/**
	 * Stop painting the model configuration on screen
	 */
	public void stopDrawingImage(){
		timer.cancel();
		timer = null;
	}
	
	@Override
	public void paint(Graphics g){
		super.paint(g);
		g.drawImage(
				fg, 0, this.getInsets().top, 
				this.getWidth(), 
				this.getHeight() - this.getInsets().top, null);
	}
	
	@Override
	public void update(int time) {
		drawImage();
	}	
	
	public void drawImage(){
	  int width = model.getLx();
    int height = model.getLy();
	  CellGroup group1 = model.getCellGroup(0);
//	  CellGroup group2 = model.getCellGroup(1);
    double value1, value2;
    int binIndex1, binIndex2;
    //draw an initial image of the cell
    for (int i = 0; i < width; i++){
      for (int j = 0; j < height; j++){
        value1 = group1.get(i, j);
//        value2 = group2.get(i, j);
        if (value1 >= maxVal) value1 = maxVal-0.00001;
        if (value1 < minVal) value1 = minVal;
        binIndex1 = (int) (value1 / incVal);
/*        if (value2 >= maxVal) value2 = maxVal-0.00001;
        if (value2 < minVal) value2 = minVal;
        binIndex2 = (int) (value2 / incVal);
        if (value1 > value2){
          fg.setRGB(i, j, colours[0][binIndex1].getRGB());      
        } else {
          fg.setRGB(i, j, colours[1][binIndex2].getRGB());
        }*/
        fg.setRGB(i, j, colours[0][binIndex1].getRGB());
      }
    }
	}
}
