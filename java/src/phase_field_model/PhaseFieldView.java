package phase_field_model;

import java.awt.*;

import javax.swing.*;

/**
 * PottsView.java
 * 
 * Handle the window operations for the GUI. Main class of the program
 * for simulations with visualisation.
 * 
 * @author Michael Chiang
 *
 */

@SuppressWarnings("serial")
public class PhaseFieldView extends JFrame {
	private PhaseFieldModel model = new PhaseFieldModel(0,0,0);
	private PhaseFieldViewPanel viewPanel;
	private PhaseFieldControlPanel controlPanel;
	
	//constructor
	public PhaseFieldView(){
		this.setSize(1000, 1000);
		this.setTitle("Phase Field Model");
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		this.setResizable(true);
		
		viewPanel = new PhaseFieldViewPanel(model);
		controlPanel = new PhaseFieldControlPanel(this);
		
		Container content = this.getContentPane();
		content.add(viewPanel, BorderLayout.CENTER);
		content.add(controlPanel, BorderLayout.SOUTH);
		
		this.setVisible(true);
	}
	
	public static void main (String [] args){
		new PhaseFieldView();
	}
	
	/**
	 * Set the model to be display on screen
	 * @param model
	 */
	public void setModel(PhaseFieldModel model){
		viewPanel.setModel(model);
	}
	
	/**
	 * Start drawing the model configuration to screen
	 */
	public void initImage(){
		viewPanel.initImage();
	}
	
	/**
	 * Stop drawing the model configuration to screen
	 */
	public void stopDrawingImage(){
		viewPanel.stopDrawingImage();
	}
}
