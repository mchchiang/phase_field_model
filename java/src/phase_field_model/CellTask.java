package phase_field_model;

import java.util.concurrent.Callable;

public class CellTask implements Callable<Object> {
  
  private PhaseFieldModel model;
  private Cell cell;
  
  public CellTask (PhaseFieldModel model, Cell cell){
    this.model = model;
    this.cell = cell;
  }
  
  @Override
  public Object call() throws Exception {
    model.updateCellVolume(cell);
    model.updateCellField(cell);
    model.updateCellCM(cell);
    return null;
  }

}
