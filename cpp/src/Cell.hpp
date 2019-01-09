/*
 * Cell.hpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#ifndef CELL_HPP_
#define CELL_HPP_

#include <vector>
#include <random>

#include "Field2D.hpp"
#include "VolumeField.hpp"

class Cell : public Field2D {

private:
  int x {}, y {};
  double px {}, py {};
  double theta {}, rotateDiff {0.01};
  int lx {}, ly {};
  double volume {};
  double xcm {}, ycm {};
  double deltaXCM {}, deltaYCM {};
  int cellType {};
  int setField {}, getField {};
  const double INCELL {0.4};
  const double CMSHIFT {2.0};

  // A random generator
  std::mt19937 mt;
  std::uniform_real_distribution<double> randDouble;

  VolumeField volumeField {VolumeField(this)};
  std::vector<std::vector<double> > cellField[2];

protected:
  double calculateTotalVolume();
  std::vector<double> calculateCM();


public:

  Cell(int x, int y, int lx, int ly, int type);

  void initOnes(int x, int y, int dx, int dy);
  void initCell(const std::vector<std::vector<double> >& matrix);
  void updateTotalVolume();
  void updateCM();
  void updateVelocity();
  void startUpdateCellField();
  void endUpdateCellField();
  void shiftCoordinates(int xShift, int yShift);

  // Accessor methods
  double getXCM() const;
  double getYCM() const;
  int getX() const;
  int getY() const;
  void setTheta(double angle);
  double getTheta() const;
  double getPx() const;
  double getPy() const;
  double getVolume(int i, int j) const;
  double getTotalVolume() const;
  Field2D* getVolumeField();
  int getCellType() const;

  //  Methods from Field2D
  int getLx() const;
  int getLy() const;
  void set(int i, int j, double value);
  double get(int i, int j) const;
};

#endif /* CELL_HPP_ */
