package phase_field_model;

public class Cell implements Field2D {

  private int lx, ly;
  private double volume;

  // Store cm coordinates in cell's own reference frame
  private double xcm, ycm;
  private double deltaXCM, deltaYCM;

  private int cellType;
  private int x, y;
  private double cellField[][][];
  private int setField, getField; // For switching between old and new field
  private static double INCELL = 0.4;
  private static double CMSHIFT = 2.0;
  private VolumeField volumeField;

  public Cell(int x, int y, int lx, int ly, int type) {
    this.x = x;
    this.y = y;
    this.lx = lx;
    this.ly = ly;
    this.xcm = lx / 2.0;
    this.ycm = ly / 2.0;
    this.deltaXCM = 0.0;
    this.deltaYCM = 0.0;
    this.cellType = type;
    this.cellField = new double[2][lx][ly]; // Default field value is zero
    this.setField = 1;
    this.getField = 1;
    this.volumeField = new VolumeField(this);
  }

  public void initSquareCell(int dx, int dy){
    int xStart = (int) (lx/2.0 - dx/2.0);
    int xEnd = xStart + dx;
    int yStart = (int) (ly/2.0 - dy/2.0);
    int yEnd = yStart + dy;
    for (int k = 0; k < 2; k++){
      for (int i = xStart; i < xEnd; i++){
        for (int j = yStart; j < yEnd; j++){
          cellField[k][i][j] = 1.0;
        }
      }
    }
    double [] cm = calculateCM();
    xcm = cm[0];
    ycm = cm[1];
  }
  
  public void initCell(double [][] matrix){
    if (matrix.length == lx && matrix[0].length == ly){
      for (int i = 0; i < lx; i++){
        for (int j = 0; j < ly; j++){
          cellField[0][i][j] = matrix[i][j];
          cellField[1][i][j] = matrix[i][j];
        }
      }
    }
    double [] cm = calculateCM();
    xcm = cm[0];
    ycm = cm[1];
  }

  public void updateTotalVolume() {
    volume = calculateTotalVolume();
  }
  
  protected double calculateTotalVolume() {
    double sum = 0.0;
    for (int i = 0; i < lx; i++) {
      for (int j = 0; j < ly; j++) {
        sum += volumeField.get(i, j);
      }
    }
    return sum;
  }
  
  protected double [] calculateCM() {
    double xavg = 0;
    double yavg = 0;
    int count = 0;
    for (int i = 0; i < lx; i++) {
      for (int j = 0; j < ly; j++) {
        if (cellField[getField][i][j] > INCELL) {
          xavg += i;
          yavg += j;
          count++;
        }
      }
    }
    if (count > 0){
      xavg /= (double) count;
      yavg /= (double) count;
    } else {
      xavg = 0.0;
      yavg = 0.0;
    }
    return new double [] {xavg, yavg};
  }

  public void updateCM() {
	  
    double oldXCM = xcm;
    double oldYCM = ycm;
    
    double [] cm = calculateCM();
    xcm = cm[0];
    ycm = cm[1];

    deltaXCM += (xcm - oldXCM);
    deltaYCM += (ycm - oldYCM);
    
    if (Math.abs(deltaXCM) > CMSHIFT || Math.abs(deltaYCM) > CMSHIFT) {
      int xShift = (int) Math.floor(deltaXCM);
      int yShift = (int) Math.floor(deltaYCM);
      shiftCoordinates(xShift, yShift);
      deltaXCM = 0;
      deltaYCM = 0;
      cm = calculateCM();
      xcm = cm[0];
      ycm = cm[1];
    }
  }

  public void shiftCoordinates(int xShift, int yShift) {
    startUpdateCellField();
    int xStart, xEnd, yStart, yEnd;
    int zeroXStart, zeroXEnd, zeroYStart, zeroYEnd;
    if (xShift >= 0) {
      xStart = xShift;
      xEnd = lx;
      zeroXStart = lx - xShift;
      zeroXEnd = lx;
    } else {
      xStart = 0;
      xEnd = lx + xShift;
      zeroXStart = 0;
      zeroXEnd = -xShift;
    }
    if (yShift >= 0) {
      yStart = yShift;
      yEnd = ly;
      zeroYStart = ly - yShift;
      zeroYEnd = ly;
    } else {
      yStart = 0;
      yEnd = ly + yShift;
      zeroYStart = 0;
      zeroYEnd = -yShift;
    }
    
    // Set empty cells to zero
    for (int i = zeroXStart; i < zeroXEnd; i++) {
      for (int j = 0; j < ly; j++) {
        set(i, j, 0.0);
      }
    }
    for (int i = zeroXEnd; i < lx; i++) {
      for (int j = zeroYStart; j < zeroYEnd; j++) {
        set(i, j, 0.0);
      }
    }
    
    // Shift cells
    for (int i = xStart; i < xEnd; i++) {
      for (int j = yStart; j < yEnd; j++) {
        set(i - xShift, j - yShift, get(i, j));
      }
    }
    
    // Update end coordinate
    x += xShift;  
    y += yShift;
    endUpdateCellField();
  }

  // Accessor methods
  public double getXCM() {
    return xcm;
  }

  public double getYCM() {
    return ycm;
  }

  public int getX() {
    return x;
  }

  public int getY() {
    return y;
  }
  
  public double getVolume(int i, int j) {
    return volumeField.get(i, j);
  }

  public double getTotalVolume() {
    return volume;
  }

  public Field2D getVolumeField() {
    return volumeField;
  }

  public int getCellType() {
    return cellType;
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
    return cellField[getField][i][j];
  }

  @Override
  public void set(int i, int j, double value) {
    cellField[setField][i][j] = value;
  }
  
  protected void startUpdateCellField() {
    setField = getField == 1 ? 0 : 1;
  }

  protected void endUpdateCellField() {
    getField = setField;
  }
}
