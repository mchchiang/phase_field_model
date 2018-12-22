/*
 * Field2D.hpp
 *
 *  Created on: 21 Dec, 2018
 *      Author: MichaelChiang
 */

#ifndef FIELD2D_HPP_
#define FIELD2D_HPP_

class Field2D {

public:
  virtual ~Field2D(){}
  virtual int getLx() const = 0;
  virtual int getLy() const = 0;
  virtual void set(int i, int j, double value) = 0;
  virtual double get(int i, int j) const = 0;

};

#endif /* FIELD2D_HPP_ */
