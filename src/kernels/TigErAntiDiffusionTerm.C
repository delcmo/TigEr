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

#include "TigErAntiDiffusionTerm.h"
/**
This function computes the artificial dissipative terms that when used with mass lumping, verifies the maximum principle. The implementation of this function is based upon the paper by Guermond and Popov: 'A Maximum-Principle Preserving C0 Finite Element Method for Scalar Conservation Equations'.
**/
template<>
InputParameters validParams<TigErAntiDiffusionTerm>()
{
  InputParameters params = validParams<Kernel>();

  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");
  // Angular:
  params.addParam<Real>("angular_direction", 1., "angular direction of the radiation: +1 or -1");

  return params;
}

TigErAntiDiffusionTerm::TigErAntiDiffusionTerm(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // get the nodal values
    _u_nodal_old(_var.nodalValueOld()),
    _u_nodal(_var.nodalValue()),
    // Speed of light constant:
    _c(getParam<Real>("c")),
    // Angular
    _omega(getParam<Real>("angular_direction")),
    // Material properties:
    _kappa(getMaterialProperty<Real>("kappa")),
    _sigma(getMaterialProperty<Real>("sigma"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

void
TigErAntiDiffusionTerm::computeResidual()
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

Real TigErAntiDiffusionTerm::computeQpResidual()
{
  // Get the nodal values and the shape functions
  VariablePhiValue phi = _var.phi();
  VariablePhiGradient phi_grad = _var.gradPhi();

  // Compute the local mass matrix M_{i,j}, A_{i,j} and b_k(\phi_j, \phi_i) in the same loop
  std::vector<std::vector<Real> > Mij(phi.size(), std::vector<Real>(phi.size(), 0.));
  std::vector<std::vector<Real> > Aij(phi.size(), std::vector<Real>(phi.size(), 0.));
  std::vector<std::vector<Real> > Dij(phi.size(), std::vector<Real>(phi.size(), 0.));
  std::vector<std::vector<Real> > b_k(phi.size(), std::vector<Real>(phi.size(), 0.));  
  std::vector<Real> b_i(phi.size(), 0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
  {
    for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
      b_i[ivar] += 0.;//_c*_test[ivar][qp];

    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      b_k[ivar][jvar] = jvar == ivar ? _current_elem_volume : -_current_elem_volume;
      Dij[ivar][jvar] = b_k[ivar][jvar];
      for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
      {
        Mij[ivar][jvar] += phi[jvar][_qp]*_test[ivar][qp]*_coord[qp]*_JxW[qp];
        Aij[ivar][jvar] = (_omega*phi_grad[jvar][qp](0)+_sigma[_qp]*phi[jvar][qp])*_c*_test[ivar][qp]*_coord[qp]*_JxW[qp];
      }
    }
  }

  // Compute the low-order viscosity coefficient
  Real sum_b_k(0.), max_Aij(0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      if (ivar != jvar)
        max_Aij = std::max(max_Aij, std::max(0., Aij[ivar][jvar])/(-b_k[ivar][jvar]));

  Real low_visc = max_Aij;

  if (low_visc<0)
    mooseError("The low-order viscosity coefficient computed in '"<<this->name()<<"' is locally negative.");

  // Compute the upper and lower bounds for the fluxes at node i: W_plus and W_minus
  Real U_minus(0.), U_plus(0.);
  for (unsigned int jvar=0; jvar<phi.size(); jvar++)
  {
    U_minus = std::min(U_minus, _u_nodal_old[jvar]);
    U_plus = std::max(U_plus, _u_nodal_old[jvar]);
  }

  Real W_plus(0.), W_minus(0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
  {
    Real W(0.);
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      W -= _dt/Mij[ivar][ivar]*(Aij[ivar][jvar]+low_visc*Dij[ivar][jvar]);

    W += 1.;
    W_minus = U_minus*W + _dt/Mij[ivar][ivar]*b_i[ivar];
    W_plus = U_plus*W + _dt/Mij[ivar][ivar]*b_i[ivar];
  }

  // Compute the antidiffusion term:
  std::vector<std::vector<Real> > Fij(phi.size(), std::vector<Real>(phi.size(), 0.));
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      Fij[ivar][jvar] = -Mij[ivar][jvar]/_dt;
      Fij[ivar][jvar] *= _u_nodal[jvar]-_u_nodal_old[jvar]-_u_nodal[ivar]+_u_nodal_old[ivar];
      Fij[ivar][jvar] += Dij[ivar][jvar]*(low_visc-_kappa[_qp])*(_u_nodal_old[jvar]-_u_nodal_old[ivar]);
    }

  // Compute the limiter coefficient alpha_ij;
  std::vector<std::vector<Real> > alpha_ij(phi.size(), std::vector<Real>(phi.size(), 0.));

  // Compute the diffusion terms Dl
  Real Dl(0.), Dh(0.);
  for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    Dl += _u_nodal[jvar]*b_k[_i][jvar];
  Dl *= low_visc;

  // Return value
  return Dl;
}

Real TigErAntiDiffusionTerm::computeQpJacobian()
{
  return 0.;
}

Real TigErAntiDiffusionTerm::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}