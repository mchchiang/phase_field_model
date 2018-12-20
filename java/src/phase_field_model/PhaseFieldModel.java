package phase_field_model;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class PhaseFieldModel {

  // Lattice spacing dx = 1
  private int lx, ly;
  private double dt = 0.005;

  private ArrayList<Double> idealCellVolume; // Ideal cell volume
  private ArrayList<Double> alpha; // Growth coefficient
  private ArrayList<Double> beta; // Excluded volume coefficient
  private ArrayList<Double> gamma; // Regularisation for adhesion
  private ArrayList<Double> eta; // Adhesion coefficient
  private ArrayList<Double> D; // Surface-tension-like coefficient
  private ArrayList<Integer> cellLx;
  private ArrayList<Integer> cellLy;

  private ArrayList<CellGroup> cellGroups;
  
  private ArrayList<PhaseFieldDataListener> listeners;
  
  private boolean running = false;
  private boolean paused = false;

  public PhaseFieldModel(int lx, int ly, int numOfCellGroups) {
    this.lx = lx;
    this.ly = ly;
    
    int numOfCellGroupsSq = numOfCellGroups * numOfCellGroups;
    cellGroups = new ArrayList<CellGroup>(numOfCellGroups); 
    for (int i = 0; i < numOfCellGroups; i++){
      cellGroups.add(new CellGroup(lx, ly));
    }
    D = new ArrayList<Double>(Collections.nCopies(numOfCellGroups, 0.0));
    alpha = new ArrayList<Double>(Collections.nCopies(numOfCellGroups, 0.0));
    idealCellVolume = 
        new ArrayList<Double>(Collections.nCopies(numOfCellGroups, 0.0));
    gamma = new ArrayList<Double>(Collections.nCopies(numOfCellGroups, 0.0));
    beta = new ArrayList<Double>(Collections.nCopies(numOfCellGroupsSq, 0.0));
    eta = new ArrayList<Double>(Collections.nCopies(numOfCellGroupsSq, 0.0)); 
    cellLx = new ArrayList<Integer>(Collections.nCopies(numOfCellGroups, 0));
    cellLy = new ArrayList<Integer>(Collections.nCopies(numOfCellGroups, 0));
    listeners = new ArrayList<PhaseFieldDataListener>();
  }
  
  public void initSquareCellLattice(int x0, int y0, int xlen, int ylen,
      int cx, int cy, int numOfCells, int type){
    double volumePerCell = xlen*ylen/(double)numOfCells;
    int dx = (int) Math.floor(Math.sqrt(volumePerCell));
    int dy = dx;
    int x = x0;
    int y = y0;
    CellGroup group = cellGroups.get(type);
    cellLx.set(type, cx);
    cellLy.set(type, cy);
    double idealVolume = idealCellVolume.get(type);
    int cellLen = (int) Math.floor(Math.sqrt(idealVolume));
    for (int i = 0; i < numOfCells; i++){
      Cell cell = new Cell(x, y, cellLx.get(type), cellLy.get(type), type);
      cell.initSquareCell(cellLen, cellLen);
      group.addCell(cell);
      x += dx;
      if (x > x0 + xlen){
        x = x0;
        y += dy;
      }
    }
  }

  public void addDataListener(PhaseFieldDataListener l){
    listeners.add(l);
  }
  
  public void removeDataListener(PhaseFieldDataListener l){
    listeners.remove(l);
  }
  
  public void notifyDataListeners(int time){
    for (PhaseFieldDataListener listener : listeners){
      listener.update(time);
    }
  }
  
  public void run(int nsteps, int threads) {
    
    running = true;
    
    ExecutorService es = Executors.newFixedThreadPool(threads);
        
    for (int i = 0; i < nsteps && running; i++) {
      updateCellGroupVolume();
//      updateCellVolume();
//      nextStep();
//      updateCellCM();
      
      List<Callable<Object>> todos = new ArrayList<Callable<Object>>();
      
      for (CellGroup group : cellGroups){
        Iterator<Cell> it = group.getCells();
        while (it.hasNext()){
          todos.add(new CellTask(this, it.next()));
        }
      }
      
      try {
        es.invokeAll(todos);
      } catch (InterruptedException e) {
        e.printStackTrace();
      }
      
      output(i);
      
    //for pausing the simulation
      synchronized(this){
        if (paused){
          try {
            this.wait();
          } catch (InterruptedException e) {}
        }
      }
    }
    
    es.shutdown();
    
    running = false;
  }

  public void updateCellGroupVolume() {
    for (CellGroup group : cellGroups) {
      group.updateField();
    }
  }

  public void updateCellVolume() {
    for (CellGroup group : cellGroups) {
      Iterator<Cell> it = group.getCells();
      while (it.hasNext()) {
        updateCellVolume(it.next());
      }
    }
  }
  
  public void updateCellVolume(Cell cell){
    cell.updateTotalVolume();
  }

  public void nextStep() {
    for (CellGroup group : cellGroups) {
      Iterator<Cell> it = group.getCells();
      while (it.hasNext()) {
        Cell cell = it.next();
        updateCellField(cell);
      }
    }
  }

  public void updateCellField(Cell cell) {
    // Apply fixed (Dirichlet) boundary conditions (u = 0 at boundaries)
    // i and j are coordinates in the cell's reference frame
    cell.startUpdateCellField();
    int type = cell.getCellType();
    int cx = cellLx.get(type);
    int cy = cellLy.get(type);
    for (int i = 1; i < cx - 1; i++) {
      for (int j = 1; j < cy - 1; j++) {
        cell.set(i, j, cell.get(i, j) + dt * (
            singleCellInteractions(cell, i, j) +
            cellCellInteractions(cell, i, j) +
            cellSubstrateInteractions(cell, i, j)));
      }
    }
    cell.endUpdateCellField();
  }

  public double singleCellInteractions(Cell cell, int i, int j) {
    double u = cell.get(i, j);
    int type = cell.getCellType();
    return D.get(type) * centralDiff(i, j, cell) + 
        0.05*(forwardDiff(i, j, 0, cell) + forwardDiff(i, j, 1, cell)) +
        u * (1 - u) * (u - 0.5 + getVolumeCoeff(type)
            * (getIdealCellVolume(type) - cell.getTotalVolume()));
  }
  
  public double cellCellInteractions(Cell cell, int i, int j) {
    double u = cell.get(i, j);
    int type = cell.getCellType();
    int cellX = cell.getX();
    int cellY = cell.getY();
    int x = iwrap(cellX + i);
    int y = jwrap(cellY + j);

    CellGroup group;

    // Compute excluded volume effect
    double exclusion = getExclusionCoeff(type, type) * cell.getVolume(i, j);
    for (int l = 0; l < cellGroups.size(); l++) {
      group = cellGroups.get(l);
      exclusion -= getExclusionCoeff(l, type) * group.get(x, y);
    }

    // Compute adhesion effect
    double volumeFieldCentralDiff = centralDiff(i, j, cell.getVolumeField());
    double adhesion = -getAdhesionCoeff(type, type) * volumeFieldCentralDiff;
    for (int l = 0; l < cellGroups.size(); l++) {
      group = cellGroups.get(l);
      adhesion += getAdhesionCoeff(l, type) * centralDiff(x, y, group);
    }

    // Compute regularisation effect
    double regularisation = getRegulateCoeff(type) * volumeFieldCentralDiff;

    return u * (1 - u) * (exclusion + adhesion + regularisation);
  }

  public double cellSubstrateInteractions(Cell cell, int i, int j) {
    return 0.0; // Not considered at the moment
  }

  public void updateCellCM() {
    for (CellGroup group : cellGroups) {
      Iterator<Cell> it = group.getCells();
      while (it.hasNext()) {
        updateCellCM(it.next());       
      }
    }
  }
  
  public void updateCellCM(Cell cell){
    cell.updateCM();
  }
  
  public void output(int step) {
    if (step % 1000 == 0) {
      notifyDataListeners(step);
      System.out.println("Step " + step);
      try {
        PrintWriter writer = new PrintWriter(new FileWriter("output.dat"));
        writer.print(cellGroups.get(0).toString());
        writer.close();
        PrintWriter cmwriter = new PrintWriter(new FileWriter("cm.dat"));
        for (CellGroup group : cellGroups){
          Iterator<Cell> it = group.getCells();
          while (it.hasNext()){
            Cell cell = it.next();
            double x = (cell.getX() + cell.getXCM()) % lx;
            double y = (cell.getY() + cell.getYCM()) % ly;
            cmwriter.println(String.format("%.5f %.5f", x, y));
          }
        }
        cmwriter.close();
      } catch (IOException e) {} // Do nothing for now
    }  
  }
  
  // For animations
  public synchronized void stop(){
    if (paused){
      resume();
    }
    running = false;
  }

  public void pause(){
    paused = true;
  }

  public synchronized void resume(){
    paused = false;
    this.notifyAll();
  }

  public boolean isRunning(){
    return running;
  }

  public boolean isPaused(){
    return paused;
  }


  // Accessor methods
  public int getLx() {
    return lx;
  }

  public int getLy() {
    return ly;
  }
  
  public int getCellLx(int type){
    return cellLx.get(type);
  }
  
  public int getCellLy(int type){
    return cellLy.get(type);
  }

  public int getNumOfCells() {
    int cells = 0;
    for (CellGroup group : cellGroups) {
      cells += group.getNumOfCells();
    }
    return cells;
  }

  public int getNumOfCellGroups() {
    return cellGroups.size();
  }
  
  public CellGroup getCellGroup(int type){
    return cellGroups.get(type);
  }
  
  public void setIdealCellVolume(int type, double value){
    if (value >= 0.0){
      idealCellVolume.set(type, value);
    }
  }

  public double getIdealCellVolume(int type) {
    return idealCellVolume.get(type);
  }
  
  public void setVolumeCoeff(int type, double value){
    if (value >= 0.0){
      alpha.set(type, value);
    }
  }

  public double getVolumeCoeff(int type) {
    return alpha.get(type);
  }
  
  public void setExclusionCoeff(int type1, int type2, double value){
    if (value >= 0.0){
      beta.set(getTypeIndex(type1, type2), value);
      beta.set(getTypeIndex(type2, type1), value);
    }
  }

  public double getExclusionCoeff(int type1, int type2) {
    return beta.get(getTypeIndex(type1, type2));
  }

  public void setAdhesionCoeff(int type1, int type2, double value){
    if (value >= 0.0){
      eta.set(getTypeIndex(type1, type2), value);
      eta.set(getTypeIndex(type2, type1), value);
    }
  }
  
  public double getAdhesionCoeff(int type1, int type2) {
    return eta.get(getTypeIndex(type1, type2));
  }
  
  public void setRegulateCoeff(int type, double value) {
    if (value >= 0.0){
      gamma.set(type, value);
    }
  }

  public double getRegulateCoeff(int type) {
    return gamma.get(type);
  }
  
  public void setDiffusionCoeff(int type, double value){
    if (value >= 0.0){
      D.set(type, value);
    }
  }
    
  public double getDiffusionCoeff(int type) {
    return D.get(type);
  }
  
  protected double forwardDiff(int i, int j, int comp, Field2D field){
    if (comp == 0){
      return 0.5*(field.get(iup(i), j) - field.get(idown(i), j));
    } 
    return 0.5*(field.get(i, jup(j)) - field.get(i, jdown(j)));
  }
  
  protected double centralDiff(int i, int j, Field2D field) {
    return field.get(iup(i), j) + field.get(idown(i), j) + field.get(i, jup(j))
        + field.get(i, jdown(j)) - 4.0 * (field.get(i, j));
  }

  // Implement periodic boundary conditions
  protected int iwrap(int i) {
    int remainder = i % lx;
    if (remainder >= 0)
      return remainder;
    return lx + remainder;
  }

  protected int jwrap(int j) {
    int remainder = j % ly;
    if (remainder >= 0)
      return remainder;
    return ly + remainder;
  }

  protected int iup(int i) {
    if (i + 1 >= lx)
      return 0;
    return i + 1;
  }

  protected int idown(int i) {
    if (i - 1 < 0)
      return lx - 1;
    return i - 1;
  }

  protected int jup(int j) {
    if (j + 1 >= ly)
      return 0;
    return j + 1;
  }

  protected int jdown(int j) {
    if (j - 1 < 0)
      return ly - 1;
    return j - 1;
  }
  
  protected int getTypeIndex(int type1, int type2){
    return type1*getNumOfCellGroups()+type2;
  }
  
  public static void main (String [] args){
    PhaseFieldModel model = new PhaseFieldModel(100, 100, 1);
    model.setDiffusionCoeff(0, 1.0);
    model.setIdealCellVolume(0, 25);
    model.setRegulateCoeff(0, 1.0);
    model.setVolumeCoeff(0, 1.0);
    model.setExclusionCoeff(0, 0, 1.0);
    model.setAdhesionCoeff(0, 0, 0.3);
    model.initSquareCellLattice(10, 10, 70, 70, 20, 20, 80, 0);
    model.run(100000, 8);
  }

}