package phase_field_model;

import java.util.ArrayList;
import java.util.Iterator;

public class CellGroup implements Field2D {

  private int lx, ly;
  private double field[][];
  private ArrayList<Cell> cells;

  public CellGroup(int lx, int ly) {
    this.lx = lx;
    this.ly = ly;
    this.field = new double[lx][ly];
    cells = new ArrayList<Cell>();
  }
  
  public void addCell(Cell cell){
    cells.add(cell);
  }

  public void updateField() {
    // reset the field to zero
    for (int i = 0; i < lx; i++) {
      for (int j = 0; j < ly; j++) {
        field[i][j] = 0.0;
      }
    }

    for (Cell cell : cells) {
      int cellLx = cell.getLx();
      int cellLy = cell.getLy();
      int cellX = cell.getX();
      int cellY = cell.getY();
      int x, y;
      for (int i = 0; i < cellLx; i++) {
        for (int j = 0; j < cellLy; j++) {
          x = iwrap(cellX + i);
          y = jwrap(cellY + j);
          field[x][y] += cell.getVolume(i, j);
        }
      }
    }
  }

  // Accessor methods
  public int getNumOfCells() {
    return cells.size();
  }

  public Iterator<Cell> getCells() {
    return cells.iterator();
  }

  @Override
  public int getLx() {
    return lx;
  }

  @Override
  public int getLy() {
    return ly;
  }

  @Override
  public double get(int i, int j) {
    return field[i][j];
  }

  @Override
  public void set(int i, int j, double value) {
    field[i][j] = value;
  }
  
  public void printField(){
    for (int i = 0; i < lx; i++){
      for (int j = 0; j < ly; j++){
        System.out.print(field[i][j] + " ");
      }
      System.out.println();
    }
  }
  
  public String toString(){
    String matrix = "";
    for (int i = 0; i < lx; i++){
      for (int j = 0; j < ly; j++){
        matrix += String.format("%.5f ", field[i][j]);
      }
      matrix += "\n";
    }
    return matrix;  
  }

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
}
