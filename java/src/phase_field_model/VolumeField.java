package phase_field_model;

public class VolumeField implements Field2D {

  private Field2D field;

  public VolumeField(Field2D field) {
    this.field = field;
  }

  @Override
  public int getLx() {
    return field.getLx();
  }

  @Override
  public int getLy() {
    return field.getLy();
  }

  @Override
  public void set(int i, int j, double value) {
    // Cannot set values to a wrapper/derived field
  }

  @Override
  public double get(int i, int j) {
    double value = field.get(i, j);
    return value * value * (3 - 2 * value);
  }

}
