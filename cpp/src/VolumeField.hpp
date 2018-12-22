/*
 * VolumeField.hpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#ifndef VOLUMEFIELD_HPP_
#define VOLUMEFIELD_HPP_

#include "Field2D.hpp"

class VolumeField : public Field2D {

private:
  Field2D* field;

public:
  VolumeField(Field2D* f) {
    field = f;
  }

  int getLx() const {
    return field->getLx();
  }

  int getLy() const {
    return field->getLy();
  }

  void set(int i, int j, double value) {
    // Do nothing
  }

  double get(int i, int j) const {
    double value = field->get(i, j);
    return value * value * (3 - 2 * value);
  }
};

#endif /* VOLUMEFIELD_HPP_ */
