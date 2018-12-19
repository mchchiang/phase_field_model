package phase_field_model;

import org.junit.Assert;
import org.junit.Test;

public class PhaseFieldModelTest {
  
  private class TestField implements Field2D {
    
    private double [][] field;
    
    public TestField(double [][] field){
      this.field = field;
    }

    @Override
    public int getLx() {
       return field.length;
    }

    @Override
    public int getLy() {
      return field[0].length;
    }

    @Override
    public void set(int i, int j, double value) {
      field[i][j] = value;
    }

    @Override
    public double get(int i, int j) {
      return field[i][j];
    }    
  }
  
  private double tol = 0.00001;

  @Test
  public void testIWrap1() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(3, model.iwrap(3));
  }
  
  @Test
  public void testIWrap2() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(8, model.iwrap(28));
  }
  
  @Test
  public void testIWrap3() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(5, model.iwrap(-15));
  }
  
  @Test
  public void testIWrap4() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(0, model.iwrap(-10));
  }
  
  @Test
  public void testJWrap1() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(6, model.jwrap(6));
  }
  
  @Test
  public void testJWrap2() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(4, model.jwrap(36));
  }
  
  @Test
  public void testJWrap3() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(1, model.jwrap(-15));
  }
  
  @Test
  public void testJWrap4() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(0, model.jwrap(-40));
  }
  
  @Test
  public void testIUp1() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(3, model.iup(2));
  }
  
  @Test
  public void testIUp2() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(0, model.iup(9));
  }
  
  @Test
  public void testIDown1() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(5, model.idown(6));
  }
  
  @Test
  public void testIDown2() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(9, model.idown(0));
  }
  
  @Test
  public void testJUp1() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(7, model.jup(6));
  }
  
  @Test
  public void testJUp2() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(0, model.jup(8));
  }
  
  @Test
  public void testJDown1() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(3, model.jdown(4));
  }
  
  @Test
  public void testJDown2() {
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    Assert.assertEquals(7, model.jdown(0));
  }

  @Test
  public void testCentralDiff1(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.centralDiff(1, 2, field);
    Assert.assertEquals(-0.6, result, tol);
  }
  
  @Test
  public void testCentralDiff2(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.centralDiff(2, 1, field);
    Assert.assertEquals(0.8, result, tol);
  }
  
  @Test
  public void testCentralDiff3(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.centralDiff(1, 1, field);
    Assert.assertEquals(1.6, result, tol);
  }
  
  @Test
  public void testCentralDiff4(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.centralDiff(2, 2, field);
    Assert.assertEquals(-1.6, result, tol);
  }
  
  @Test
  public void testFowardDiff1X(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(1, 1, 0, field);
    Assert.assertEquals(-0.15, result, tol);
  }
  
  @Test
  public void testFowardDiff1Y(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(1, 1, 1, field);
    Assert.assertEquals(0.05, result, tol);
  }
  
  @Test
  public void testFowardDiff2X(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(1, 2, 0, field);
    Assert.assertEquals(0.25, result, tol);
  }
  
  @Test
  public void testFowardDiff2Y(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(1, 2, 1, field);
    Assert.assertEquals(0.35, result, tol);
  }
  
  @Test
  public void testFowardDiff3X(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(2, 1, 0, field);
    Assert.assertEquals(0.15, result, tol);
  }
  
  @Test
  public void testFowardDiff3Y(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(2, 1, 1, field);
    Assert.assertEquals(-0.05, result, tol);
  }
  
  @Test
  public void testFowardDiff4X(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(2, 2, 0, field);
    Assert.assertEquals(-0.2, result, tol);
  }
  
  @Test
  public void testFowardDiff4Y(){
    
    double [][] matrix = 
      {{0.2, 0.5, 0.1, 0.8},
       {0.4, 0.0, 0.5, 0.7},
       {0.7, 0.2, 0.6, 0.0},
       {1.0, 0.3, 0.1, 0.3}};
    
    TestField field = new TestField(matrix);
    
    PhaseFieldModel model = new PhaseFieldModel(10, 8, 1);
    double result = model.forwardDiff(2, 2, 1, field);
    Assert.assertEquals(-0.1, result, tol);
  }
  
  @Test
  public void testSetDiffusionCoeff1(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setDiffusionCoeff(1, 2.0);
    Assert.assertEquals(2.0, model.getDiffusionCoeff(1), tol);
  }
  
  @Test
  public void testSetDiffusionCoeff2(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setDiffusionCoeff(1, -1.0);
    Assert.assertEquals(0.0, model.getDiffusionCoeff(1), tol);
  }
  
  @Test
  public void testSetVolumeCoeff1(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setVolumeCoeff(1, 0.4);
    Assert.assertEquals(0.4, model.getVolumeCoeff(1), tol);
  }
  
  @Test
  public void testSetVolumeCoeff2(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setVolumeCoeff(1, -0.3);
    Assert.assertEquals(0.0, model.getVolumeCoeff(1), tol);
  }
  
  @Test
  public void testSetIdealCellVolume1(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setIdealCellVolume(1, 31.4);
    Assert.assertEquals(31.4, model.getIdealCellVolume(1), tol);
  }
  
  @Test
  public void testSetIdealCellVolume2(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setIdealCellVolume(1, -9.32);
    Assert.assertEquals(0.0, model.getIdealCellVolume(1), tol);
  }
  
  @Test
  public void testSetRegulateCoeff1(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setRegulateCoeff(1, 0.13);
    Assert.assertEquals(0.13, model.getRegulateCoeff(1), tol);
  }
  
  @Test
  public void testSetRegulateCoeff2(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 2);
    model.setRegulateCoeff(1, -2.3);
    Assert.assertEquals(0.0, model.getRegulateCoeff(1), tol);
  }
  
  @Test
  public void testSetAdhesionCoeff1(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setAdhesionCoeff(0, 0, 1.33);
    Assert.assertEquals(1.33, model.getAdhesionCoeff(0, 0), tol);
  }
  
  @Test
  public void testSetAdhesionCoeff2(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setAdhesionCoeff(0, 0, -3.9);
    Assert.assertEquals(0.0, model.getAdhesionCoeff(0, 0), tol);
  }
  
  @Test
  public void testSetAdhesionCoeff3(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setAdhesionCoeff(1, 2, 3.329);;
    Assert.assertEquals(3.329, model.getAdhesionCoeff(1, 2), tol);
  }
  
  @Test
  public void testSetAdhesionCoeff4(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setAdhesionCoeff(1, 2, 3.329);;
    Assert.assertEquals(3.329, model.getAdhesionCoeff(2, 1), tol);
  }
  
  @Test
  public void testSetExclusionCoeff1(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setExclusionCoeff(0, 0, 2.56);
    Assert.assertEquals(2.56, model.getExclusionCoeff(0, 0), tol);
  }
  
  @Test
  public void testSetExclusionCoeff2(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setExclusionCoeff(0, 0, -1.2);
    Assert.assertEquals(0.0, model.getExclusionCoeff(0, 0), tol);
  }  
  
  @Test
  public void testSetExclusionCoeff3(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setExclusionCoeff(2, 1, 6.13);
    Assert.assertEquals(6.13, model.getExclusionCoeff(2, 1), tol);
  }
  
  @Test
  public void testSetExclusionCoeff4(){
    PhaseFieldModel model = new PhaseFieldModel(3, 4, 3);
    model.setExclusionCoeff(2, 1, 6.13);
    Assert.assertEquals(6.13, model.getExclusionCoeff(1, 2), tol);
  }
  
  @Test
  public void testGetNumOfCellGroups1(){
    PhaseFieldModel model = new PhaseFieldModel(3, 2, 8);
    Assert.assertEquals(8, model.getNumOfCellGroups());
  }
  
  @Test
  public void testGetNumOfCellGroups2(){
    PhaseFieldModel model = new PhaseFieldModel(3, 2, 0);
    Assert.assertEquals(0, model.getNumOfCellGroups());
  }
  
  
}
