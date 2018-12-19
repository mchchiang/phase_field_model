package phase_field_model;

import org.junit.Assert;
import org.junit.Test;

public class CellTest {
  
  private double tol = 0.00001;
  
  @Test
  public void testGet() {
    Cell cell = new Cell(0, 0, 5, 5, 0);
    Assert.assertEquals(0.0, cell.get(2,3), tol);
  }

  @Test
  public void testSet() {
    Cell cell = new Cell(0, 0, 5, 5, 0);
    double [][] matrix = 
        {{0, 0.6, 0, 0.4, 0},
        {0, 0.2, 1, 1, 0.5},
        {0, 1, 0.3, 0.2, 0},
        {0, 0.5, 1, 0.8, 0},
        {0, 0, 0.9, 0, 1.2}};
    
    for (int i = 0; i < 5; i++){
      for (int j = 0; j < 5; j++){
        cell.set(i, j, matrix[i][j]);
      }
    }
    
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++){
        Assert.assertEquals(matrix[i][j], cell.get(i, j), tol);
      }
    }
  }
  
  @Test
  public void testInitCell(){
    Cell cell = new Cell(0, 0, 5, 5, 0);
    double [][] matrix = 
      {{0, 0.6, 0, 0.4, 0},
      {0, 0.2, 1, 1, 0.5},
      {0, 1, 0.3, 0.2, 0},
      {0, 0.5, 1, 0.8, 0},
      {0, 0, 0.9, 0, 1.2}};
    cell.initCell(matrix);
    
    for (int i = 0; i < 5; i++){
      for (int j = 0; j < 5; j++){
        Assert.assertEquals(matrix[i][j], cell.get(i, j), tol);
      }
    }
  }
  
  
  @Test
  public void testGetX(){
    Cell cell = new Cell(-4, 3, 5, 5, 0);
    Assert.assertEquals(-4, cell.getX());
  }
  
  @Test
  public void testGetY(){
    Cell cell = new Cell(-4, 3, 5, 5, 0);
    Assert.assertEquals(3, cell.getY());
  }
  
  @Test
  public void testCalculuateCM1(){
    Cell cell = new Cell(0, 0, 3, 3, 0);
    
    double [][] matrix = 
      {{0, 0, 0},
       {0, 0, 0},
       {0, 0, 0}};
    
    double [] expectedCM = {0.0, 0.0};
    
    cell.initCell(matrix);
    
    Assert.assertArrayEquals(expectedCM, cell.calculateCM(), tol);
  }
  
  @Test
  public void testCalculcateCM2(){
    Cell cell = new Cell(0, 0, 3, 3, 0);
    
    double [][] matrix = 
      {{0.2, 0.4, 0.7},
       {0.1, 1.0, 0.8},
       {0.3, 0.6, 0.0}};
    
    double [] expectedCM = {1.0, 1.5};
    
    cell.initCell(matrix);
    
    Assert.assertArrayEquals(expectedCM, cell.calculateCM(), tol);
  }
  
  @Test
  public void testGetVolume1() {
    Cell cell = new Cell(0, 0, 3, 3, 0);
    
    double [][] matrix = 
      {{0, 0, 0},
       {0, 0, 0},
       {0, 0, 0}};
    
    double [][] expected = 
      {{0, 0, 0},
       {0, 0, 0},
       {0, 0, 0}};
    
    cell.initCell(matrix);
    
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        Assert.assertEquals(expected[i][j], cell.getVolume(i, j), tol);
      }
    }
  }
  
  @Test
  public void testGetVolume2(){
    Cell cell = new Cell(0, 0, 3, 3, 0);
    
    double [][] matrix = 
      {{0.2, 0.4, 0.7},
       {0.1, 1.0, 0.8},
       {0.3, 0.6, 0.0}};
    
    double [][] expected = 
      {{0.104, 0.352, 0.784},
       {0.028, 1.000, 0.896},
       {0.216, 0.648, 0.000}};
    
    cell.initCell(matrix);
    
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        Assert.assertEquals(expected[i][j], cell.getVolume(i, j), tol);
      }
    }
  }
  
  @Test
  public void testCalculateTotalVolume1(){
    Cell cell = new Cell(0, 0, 3, 3, 0);
    
    double [][] matrix = 
      {{0, 0, 0},
       {0, 0, 0},
       {0, 0, 0}};
    
    cell.initCell(matrix);
    
    Assert.assertEquals(0.0, cell.calculateTotalVolume(), tol);
  }
  
  @Test
  public void testCalculateTotalVolume2(){
    Cell cell = new Cell(0, 0, 3, 3, 0);
    
    double [][] matrix = 
      {{0.2, 0.4, 0.7},
       {0.1, 1.0, 0.8},
       {0.3, 0.6, 0.0}};
    
    cell.initCell(matrix);
    
    Assert.assertEquals(4.028, cell.calculateTotalVolume(), tol);
  }
  
  @Test
  public void testShiftCoordinates1(){
    Cell cell = new Cell(0, 0, 5, 5, 0);
    double [][] matrix = 
        {{0, 0, 0, 0, 0},
          {0, 1, 0.5, 1, 0},
          {0, 0.5, 1, 1, 0},
          {0, 1, 1, 0.5, 0},
          {0, 0, 0, 0, 0}};
    
    double [][] expected = 
      {{1, 0.5, 1, 0, 0},
        {0.5, 1, 1, 0, 0},
        {1, 1, 0.5, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}
      };
    
    cell.initCell(matrix);
    cell.shiftCoordinates(1, 1);
    
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++){
        Assert.assertEquals(expected[i][j], cell.get(i, j), tol);
      }
    }
  }
  
  @Test
  public void testShiftCoordinates2(){
    Cell cell = new Cell(0, 0, 5, 5, 0);
    double [][] matrix = 
      {{0, 0, 0, 0, 0},
        {0.5, 1, 0.5, 1, 0},
        {1, 0.5, 1, 1, 0},
        {1, 1, 1, 0.5, 0},
        {0, 0, 0, 0, 0}};
    
    double [][] expected = 
      {{0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {1, 0.5, 1, 0, 0},
        {0.5, 1, 1, 0, 0},
        {1, 1, 0.5, 0, 0}
      };
    
    cell.initCell(matrix);   
    cell.shiftCoordinates(-1, 1);
    
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++){
        Assert.assertEquals(expected[i][j], cell.get(i, j), tol);
      }
    }
  }
}
