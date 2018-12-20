package phase_field_model;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 * PhaseFieldControlPanel.java
 * 
 * The parameter control panel on the GUI. It processes the parameter
 * inputs by the user and initialise the simulation when the "run"
 * button is clicked.
 * 
 * @author Michael Chiang
 *
 */
@SuppressWarnings("serial")
public class PhaseFieldControlPanel extends JPanel implements ActionListener {
	private JPanel simParamsPanel;
	private JPanel modelParamsPanel;
	
//	private JLabel lblDiffusionCoeff;
//	private JTextField txtDiffusionCoeff;
	
	private JLabel lblNumOfCells;
	private JTextField txtNumOfCells;
	
	private JLabel lblIdealVolume;
	private JTextField txtIdealVolume;
	
	private JLabel lblExclusionCoeff;
	private JTextField txtExclusionCoeff;
	
	private JLabel lblAdhesionCoeff;
	private JTextField txtAdhesionCoeff;
	
	private JLabel lblCellWidth;
	private JTextField txtCellWidth;
	
	private JLabel lblCellHeight;
	private JTextField txtCellHeight;
	
	private JLabel lblWidth;
	private JTextField txtWidth;
	
	private JLabel lblHeight;
	private JTextField txtHeight;
	
	private JLabel lblNumOfThreads;
	private JTextField txtNumOfThreads;
	
	private JLabel lblNumOfSteps;
	private JTextField txtNumOfSteps;
	
	private JLabel lblDt;
	private JTextField txtDt;
	
	private JButton btnRun;
	private JButton btnStop;
	private JButton btnPause;
	
	private PhaseFieldView view;
	private PhaseFieldModel model;
	
	public PhaseFieldControlPanel(PhaseFieldView view){
		this.view = view;
		
		lblWidth = new JLabel("Width: ");
		txtWidth = new JTextField(3);
		
		lblHeight = new JLabel("Height: ");
		txtHeight = new JTextField(3);
		
		lblNumOfThreads = new JLabel("NThreads: ");
		txtNumOfThreads = new JTextField(3);
		
		lblDt = new JLabel("dt: ");
    txtDt = new JTextField(3);
		
		txtNumOfSteps = new JTextField(3);
		lblNumOfSteps = new JLabel("Steps: ");
		
		btnRun = new JButton("Run");
		btnRun.addActionListener(this);
		
		btnStop = new JButton("Stop");
    btnStop.addActionListener(this);
    btnStop.setEnabled(false);

    btnPause = new JButton("Pause");
    btnPause.addActionListener(this);
    btnPause.setEnabled(false);
		
		lblNumOfCells = new JLabel("NCells: ");
		txtNumOfCells = new JTextField(3);
		
		lblCellWidth = new JLabel("Cell width: ");
		txtCellWidth = new JTextField(3);
		
		lblCellHeight = new JLabel("Cell height: ");
		txtCellHeight = new JTextField(3);
		
    lblIdealVolume = new JLabel("Vol: ");
    txtIdealVolume = new JTextField(3);   
		
		lblExclusionCoeff = new JLabel("beta: ");
		txtExclusionCoeff = new JTextField(3);
		
		lblAdhesionCoeff = new JLabel("eta: ");
		txtAdhesionCoeff = new JTextField(3);
		
		modelParamsPanel = new JPanel();
		
		modelParamsPanel.add(lblNumOfCells);
		modelParamsPanel.add(txtNumOfCells);
		modelParamsPanel.add(lblCellWidth);
		modelParamsPanel.add(txtCellWidth);
		modelParamsPanel.add(lblCellHeight);
		modelParamsPanel.add(txtCellHeight);
		modelParamsPanel.add(lblIdealVolume);
		modelParamsPanel.add(txtIdealVolume);
		modelParamsPanel.add(lblExclusionCoeff);
    modelParamsPanel.add(txtExclusionCoeff);
    modelParamsPanel.add(lblAdhesionCoeff);
    modelParamsPanel.add(txtAdhesionCoeff);
		
		simParamsPanel = new JPanel();
		simParamsPanel.add(lblWidth);
		simParamsPanel.add(txtWidth);
		simParamsPanel.add(lblHeight);
		simParamsPanel.add(txtHeight);
		simParamsPanel.add(lblDt);
		simParamsPanel.add(txtDt);
		simParamsPanel.add(lblNumOfSteps);
		simParamsPanel.add(txtNumOfSteps);
		simParamsPanel.add(lblNumOfThreads);
    simParamsPanel.add(txtNumOfThreads);
		simParamsPanel.add(btnRun);
    simParamsPanel.add(btnStop);
    simParamsPanel.add(btnPause);		
		
		setLayout(new BorderLayout());
		add(modelParamsPanel, BorderLayout.NORTH);
		add(simParamsPanel, BorderLayout.SOUTH);
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		/*
		 * run the simulation with the parameters entered by the user 
		 * when the "run" button is clicked
		 */
		if (e.getSource() == btnRun){
			Thread runthread = new Thread(){
				@Override
				public void run(){					
					int nx = Integer.parseInt(txtWidth.getText());
					int ny = Integer.parseInt(txtHeight.getText());
					int nThreads = Integer.parseInt(txtNumOfThreads.getText());
					double dt = Double.parseDouble(txtDt.getText());
					
					int numOfCells = Integer.parseInt(txtNumOfCells.getText());
					int cellLx = Integer.parseInt(txtCellWidth.getText());
					int cellLy = Integer.parseInt(txtCellHeight.getText());
					double idealVolume = Double.parseDouble(txtIdealVolume.getText());
					double exclusionCoeff = Double.parseDouble(txtExclusionCoeff.getText());
					double adhesionCoeff = Double.parseDouble(txtAdhesionCoeff.getText());
					
				  int numOfCellGroups = 1;
				  int buffer = 5;
				  
					int nsteps = Integer.parseInt(txtNumOfSteps.getText());
					
					model = new PhaseFieldModel(nx, ny, numOfCellGroups);
					model.setDiffusionCoeff(0, 1.0);
			    model.setIdealCellVolume(0, idealVolume);
			    model.setRegulateCoeff(0, 1.0);
			    model.setVolumeCoeff(0, 1.0);
			    model.setExclusionCoeff(0, 0, exclusionCoeff);
			    model.setAdhesionCoeff(0, 0, adhesionCoeff);
			    model.setNumOfThreads(nThreads);
			    model.setDt(dt);
			    model.initSquareCellLattice(10, 10, 70, 70, 
			        cellLx, cellLy, numOfCells, 0);
			    
					view.setModel(model);
					view.initImage();
					
					btnRun.setEnabled(false);
					btnStop.setEnabled(true);
					btnPause.setEnabled(true);
					
					model.run(nsteps);
					
					view.stopDrawingImage();
					
					btnRun.setEnabled(true);
				}
			};
			runthread.start();
		} else if (e.getSource() == btnStop){
      btnPause.setEnabled(false);
      btnPause.setText("Pause");
      model.stop();
      
    } else if (e.getSource() == btnPause){
      if (model.isPaused()){
        model.resume();
        btnPause.setText("Pause");
      } else {
        btnPause.setText("Resume");
        model.pause();
      }
    }
	}
}
