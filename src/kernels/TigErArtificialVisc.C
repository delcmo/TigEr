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

#include "TigErArtificialVisc.h"
/**
This function computes the artificial dissipative terms that when used with mass lumping, verifies the maximum principle. The implementation of this function is based upon the paper by Guermond and Popov: 'A Maximum-Principle Preserving C0 Finite Element Method for Scalar Conservation Equations'.
**/
template<>
InputParameters validParams<TigErArtificialVisc>()
{
  InputParameters params = validParams<Kernel>();

  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");
  // Angular:
  params.addParam<Real>("angular_direction", 1., "angular direction of the radiation: +1 or -1");

  return params;
}

TigErArtificialVisc::TigErArtificialVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // get the nodal values
    _u_nodal(_is_implicit ? _var.nodalValue() : _var.nodalValueOld()),
    // Speed of light constant:
    _c(getParam<Real>("c")),
    // Angular
    _omega(getParam<Real>("angular_direction")),
    // Material properties:
    _sigma(getMaterialProperty<Real>("sigma"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

void
TigErArtificialVisc::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); _i++)
    _local_re(_i) = computeQpResidual();

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i=0; i<_save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}

Real TigErArtificialVisc::computeQpResidual()
{
  std::cout<<_stage_fct<<std::endl;
  // Get the nodal values and the shape functions
  VariablePhiValue phi = _var.phi();
  VariablePhiGradient phi_grad = _var.gradPhi();

  // Compute A_{i,j}, component (i,j) of the steady-state matrix
  std::vector<std::vector<Real> > Aij(phi.size(), std::vector<Real>(phi.size(), 0.));
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
        Aij[ivar][jvar] += (_omega*phi_grad[jvar][qp](0)+_sigma[_qp]*phi[jvar][qp])*_c*_test[ivar][qp]*_coord[qp]*_JxW[qp];

  // Compute the bilinear form b_k(\phi_j, \phi_i)
  std::vector<std::vector<Real> > b_k(phi.size(), std::vector<Real>(phi.size(), 0.));
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      b_k[ivar][jvar] = jvar == ivar ? _current_elem_volume : -_current_elem_volume;

  // Compute the viscosity term
  Real sum_b_k(0.), max_Aij(0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      if (ivar != jvar)
        max_Aij = std::max(max_Aij, std::max(0., Aij[ivar][jvar])/(-b_k[ivar][jvar]));

  Real visc = max_Aij;
  
  if (visc<0)
    mooseError("The low-order viscosity coefficient computed in '"<<this->name()<<"' is locally negative.");

  // Compute the diffusion term
  Real diff_term(0.);
  for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    diff_term += _u_nodal[jvar]*b_k[_i][jvar];
  diff_term *= visc;

  // Return value
  return diff_term;
}

Real TigErArtificialVisc::computeQpJacobian()
{
  return 0.;
}

Real TigErArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}