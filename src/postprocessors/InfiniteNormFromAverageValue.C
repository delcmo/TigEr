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

#include "InfiniteNormFromAverageValue.h"
/* This pps computes the maximum absolute value at the quadrature points (infinite norm). */
template<>
InputParameters validParams<InfiniteNormFromAverageValue>()
{
  InputParameters params = validParams<ElementPostprocessor>();

  // Coupled variable
  params.addRequiredCoupledVar("variable", "variable this pps is acting on.");
  // Name of the postprocessor computing the average value
  params.addRequiredParam<std::string>("name_pps_average", "name of the postprocessor computing the average value.");

  return params;
}

InfiniteNormFromAverageValue::InfiniteNormFromAverageValue(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    // Coupled variable
    _u(coupledValue("variable")),
    // value of the pps
    _pps_name(getParam<std::string>("name_pps_average")),
    // Initialize the value
    _value(-std::numeric_limits<Real>::max())
{}

void
InfiniteNormFromAverageValue::initialize()
{
  _value = -std::numeric_limits<Real>::max();
  _pps_value = getPostprocessorValueByName(_pps_name);
}

void
InfiniteNormFromAverageValue::execute()
{
  // Compute maximum value within the current element:
  Real _local_max = 0.;
  for (int _qp=0; _qp < _qrule->n_points(); _qp++)
    _local_max = std::max(std::fabs(_u[_qp]-_pps_value), _local_max);

  // Get the maximum value
  _value = std::max(_value, _local_max);
}

void
InfiniteNormFromAverageValue::finalize()
{
  gatherMax(_value);
}

Real
InfiniteNormFromAverageValue::getValue()
{
  gatherMax(_value);
  return _value;
}

void
InfiniteNormFromAverageValue::threadJoin(const UserObject & y)
{
  ElementPostprocessor::threadJoin(y);
  const InfiniteNormFromAverageValue & pps = dynamic_cast<const InfiniteNormFromAverageValue &>(y);
  _value = std::max(_value, pps._value);
}
