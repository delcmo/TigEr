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

#ifndef GETEXTREMUMVALUEFROMNEIGHBORS_H
#define GETEXTREMUMVALUEFROMNEIGHBORS_H

#include "ElementPostprocessor.h"

//Forward Declarations
class GetExtremumValueFromNeighbors;

template<>
InputParameters validParams<GetExtremumValueFromNeighbors>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */

class GetExtremumValueFromNeighbors : public ElementUserObject
{
public:
  GetExtremumValueFromNeighbors(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y){}
  virtual Real getValue(){return 0.;}

  virtual void finalize();

protected:
  // Nonlinear set of variables
  NonlinearSystem & _nlsys;

  // Auxiliary set of variables
  AuxiliarySystem & _aux;

  // Variable storing the nodal extremum value
  MooseVariable & _var;
  MooseVariable & _var_out;

  // Boolean for maximum/minimum function
  bool _comp_max;
};

#endif // GETEXTREMUMVALUEFROMNEIGHBORS_H
