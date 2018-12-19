package phase_field_model;

public interface Field2D {
  public int getLx();
  public int getLy();
  public void set(int i, int j, double value);
  public double get(int i, int j);
}
