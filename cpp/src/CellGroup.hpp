/*
 * CellGroup.hpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#ifndef CELLGROUP_HPP_
#define CELLGROUP_HPP_

#include <vector>

#include "Field2D.hpp"
#include "Cell.hpp"

using std::vector;

class CellGroup : public Field2D {

private:
  int lx {}, ly {};
  std::vector<std::vector<double> > field;
  std::vector<Cell*> cells;

protected:
  int iwrap(int i);
  int jwrap(int j);

public:
  CellGroup(int lx, int ly);
  ~CellGroup();

  void addCell(Cell* cell);
  void removeCell(Cell* cell);

  void updateField();

  // Accessor methods
  int getNumOfCells() const;

  // Methods from Field2D
  int getLx() const;
  int getLy() const;
  void set(int i, int j, double value);
  double get(int i, int j) const;
};

#endif /* CELLGROUP_HPP_ */
