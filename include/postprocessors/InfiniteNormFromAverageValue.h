/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef INFINITENORMFROMAVERAGEVALUE_H
#define INFINITENORMFROMAVERAGEVALUE_H

#include "ElementPostprocessor.h"

//Forward Declarations
class InfiniteNormFromAverageValue;

template<>
InputParameters validParams<InfiniteNormFromAverageValue>();

class InfiniteNormFromAverageValue : public ElementPostprocessor
{
public:
  InfiniteNormFromAverageValue(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & y);

protected:
  // Variable this pps is acting on:
  VariableValue & _u;

  // Name of the pps
  std::string _pps_name;

  // Value of the pps
  Real _pps_value;

  // Value storing the maximum value:
  Real _value;
};

#endif // INFINITENORMFROMAVERAGEVALUE_H
