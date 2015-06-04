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

#include "ExplicitEulerFCT.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"

template<>
InputParameters validParams<ExplicitEulerFCT>()
{
  InputParameters params = validParams<TimeIntegrator>();

  return params;
}

ExplicitEulerFCT::ExplicitEulerFCT(const std::string & name, InputParameters parameters) :
    TimeIntegrator(name, parameters)
{
}

ExplicitEulerFCT::~ExplicitEulerFCT()
{
}

void
ExplicitEulerFCT::preSolve()
{
}

void
ExplicitEulerFCT::computeTimeDerivatives()
{
  _u_dot  = *_solution;
  _u_dot -= _solution_old;
  _u_dot *= 1 / _dt;
  _u_dot.close();

  _du_dot_du = 1.0 / _dt;
}

void
ExplicitEulerFCT::solve()
{
  // First stage: solve for the high-order solution
  _stage_fct = 1;
  _fe_problem.getNonlinearSystem().sys().solve();

  // Second stage: solve for the solution at next time step
}

void
ExplicitEulerFCT::postStep(NumericVector<Number> & residual)
{
  residual += _Re_time;
  residual += _Re_non_time;
  residual.close();
}
