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

#ifndef EXPLICITEULERFCT_H
#define EXPLICITEULERFCT_H

#include "TimeIntegrator.h"

class ExplicitEulerFCT;

template<>
InputParameters validParams<ExplicitEulerFCT>();

/**
 * Explicit Euler time integrator with Flux Corrected Transport (FCT)
 */
class ExplicitEulerFCT : public TimeIntegrator
{
public:
  ExplicitEulerFCT(const std::string & name, InputParameters parameters);
  virtual ~ExplicitEulerFCT();

  virtual int order() { return 1; }
  virtual void preSolve();  
  virtual void computeTimeDerivatives();
  virtual void solve();  
  virtual void postStep(NumericVector<Number> & residual);

protected:
  int _stage_fct;
};

#endif /* EXPLICITEULERFCT_H */