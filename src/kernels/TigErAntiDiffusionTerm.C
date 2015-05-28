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

  // Coupled aux variables
  params.addRequiredCoupledVar("max_nodal_values", "variable storing the maximum nodal values");
  params.addRequiredCoupledVar("min_nodal_values", "variable storing the minimum nodal values");
  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");
  // Angular:
  params.addParam<Real>("angular_direction", 1., "angular direction of the radiation: +1 or -1");

  return params;
}

TigErAntiDiffusionTerm::TigErAntiDiffusionTerm(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // get the nodal values of the variable the object is acting on
    _u_nodal_old(_var.nodalValueOld()),
    _u_nodal(_var.nodalValue()),
    // coupled aux variables: U_plus and U_minus
    _U_plus(coupledNodalValue("max_nodal_values")),
    _U_minus(coupledNodalValue("min_nodal_values")),
    // Speed of light constant:
    _c(getParam<Real>("c")),
    // Angular
    _omega(getParam<Real>("angular_direction")),
    // Material properties:
    _kappa(getMaterialProperty<Real>("entropy_viscosity_coefficient")),
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
  // Get the shape functions
  VariablePhiValue phi = _var.phi();
  VariablePhiGradient phi_grad = _var.gradPhi();

  // Compute the local mass matrix M_{i,j}, A_{i,j} Dij{i,j} and b_k(\phi_j, \phi_i) in the same loop
  std::vector<std::vector<Real> > Mij(phi.size(), std::vector<Real>(phi.size(), 0.)); // consistent mass matrix
  std::vector<std::vector<Real> > Aij(phi.size(), std::vector<Real>(phi.size(), 0.)); // stead-state transport operator
  std::vector<std::vector<Real> > Dij(phi.size(), std::vector<Real>(phi.size(), 0.)); // diffusion operator
  std::vector<std::vector<Real> > b_k(phi.size(), std::vector<Real>(phi.size(), 0.)); // bilinear form
  std::vector<Real> B_i(phi.size(), 0.); // source term vector
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
  {
    for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
      B_i[ivar] += 0.;//_c*_test[ivar][qp];

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

  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      std::cout<<"Aij["<<ivar<<"]["<<jvar<<"]="<<Aij[ivar][jvar]<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      std::cout<<"Mij["<<ivar<<"]["<<jvar<<"]="<<Mij[ivar][jvar]<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      std::cout<<"Dij["<<ivar<<"]["<<jvar<<"]="<<Dij[ivar][jvar]<<std::endl;

  // Compute the low-order viscosity coefficient
  Real sum_b_k(0.), max_Aij(0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      if (ivar != jvar)
        max_Aij = std::max(max_Aij, std::max(0., Aij[ivar][jvar])/(-b_k[ivar][jvar]));

  Real low_visc = max_Aij;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  std::cout<<"low_visc="<<low_visc<<std::endl;

  if (low_visc<0)
    mooseError("The low-order viscosity coefficient computed in '"<<this->name()<<"' is locally negative.");

  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    std::cout<<"Uminus["<<ivar<<"]="<<_U_minus[ivar]<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    std::cout<<"Uplus["<<ivar<<"]="<<_U_plus[ivar]<<std::endl;

  // Compute the upper and lower bounds for the fluxes at node i: W_plus and W_minus
  std::vector<Real> U_minus(phi.size(), 0.);
  std::vector<Real> U_plus(phi.size(), 0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
  {
    U_minus[ivar] = std::min(_U_minus[ivar], _u_nodal_old[ivar]);
    U_plus[ivar] = std::max(_U_plus[ivar], _u_nodal_old[ivar]);
  }

  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    std::cout<<"Uminus["<<ivar<<"]="<<U_minus[ivar]<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    std::cout<<"Uplus["<<ivar<<"]="<<U_plus[ivar]<<std::endl;

  std::vector<Real> W_plus(phi.size(), 0.);
  std::vector<Real> W_minus(phi.size(), 0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
  {
    Real W(0.);
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      W -= _dt/Mij[ivar][ivar]*(Aij[ivar][jvar]+low_visc*Dij[ivar][jvar]);

    W += 1.;
    W_minus[ivar] = U_minus[ivar]*W + _dt/Mij[ivar][ivar]*B_i[ivar];
    W_plus[ivar] = U_plus[ivar]*W + _dt/Mij[ivar][ivar]*B_i[ivar];
  }

  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    std::cout<<"Wminus["<<ivar<<"]="<<W_minus[ivar]<<std::endl;
  std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    std::cout<<"Wplus["<<ivar<<"]="<<W_plus[ivar]<<std::endl;

  // Compute the antidiffusion term:
  std::vector<std::vector<Real> > Fij(phi.size(), std::vector<Real>(phi.size(), 0.));
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      Fij[ivar][jvar] = -Mij[ivar][jvar]/_dt;
      Fij[ivar][jvar] *= _u_nodal[jvar]-_u_nodal_old[jvar]-_u_nodal[ivar]+_u_nodal_old[ivar];
      Fij[ivar][jvar] += Dij[ivar][jvar]*(low_visc-_kappa[_qp])*(_u_nodal_old[jvar]-_u_nodal_old[ivar]);
    }

  // Compute the sums of positive and negative fluxes: P_plus and P_minus.
  std::vector<Real> P_plus(phi.size(), 0.);
  std::vector<Real> P_minus(phi.size(), 0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      P_plus[ivar] += std::max(0., Fij[ivar][jvar]);
      P_plus[ivar] += std::min(0., Fij[ivar][jvar]);
    }

  // Compute the maximum and minimum fluxes allowed: Q_plus and Q_minus
  std::vector<Real> Q_plus(phi.size(), 0.);
  std::vector<Real> Q_minus(phi.size(), 0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
  {
    Q_plus[ivar] = Mij[ivar][ivar]*(W_plus[ivar]-_u_nodal_old[ivar])/_dt;
    Q_minus[ivar] = Mij[ivar][ivar]*(W_minus[ivar]-_u_nodal_old[ivar])/_dt;
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      Q_plus[ivar] += (Aij[ivar][jvar]+low_visc*Dij[ivar][jvar])*_u_nodal_old[jvar];
      Q_minus[ivar] += (Aij[ivar][jvar]+low_visc*Dij[ivar][jvar])*_u_nodal_old[jvar];
    }
  }

  // Compute the coefficients R_plus and R_minus
  std::vector<Real> R_plus(phi.size(), 0.);
  std::vector<Real> R_minus(phi.size(), 0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
  {
    if (P_plus[ivar] > 1.e-6)
      R_plus[ivar] = std::min(1., Q_plus[ivar]/P_plus[ivar]);

    if (P_minus[ivar] > 1.e-6)
      R_minus[ivar] = std::min(1., Q_minus[ivar]/P_minus[ivar]);
  }

  // Compute the limiter coefficient alpha_ij;
  std::vector<std::vector<Real> > alpha_ij(phi.size(), std::vector<Real>(phi.size(), 0.));
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      if (Fij[ivar][jvar] >= 0)
        alpha_ij[ivar][jvar] = std::min(R_plus[ivar], R_minus[jvar]);
      else
        alpha_ij[ivar][jvar] = std::min(R_minus[ivar], R_plus[jvar]);
    }

  // Compute the antidiffusion flux: L[F])_i = sum_j alpha_ij Fij
  std::vector<Real> antidiffusion(phi.size(), 0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      antidiffusion[ivar] += alpha_ij[ivar][jvar]*Fij[ivar][jvar];

  // Return antidiffusion flux for node i
  return antidiffusion[_i];
}

Real TigErAntiDiffusionTerm::computeQpJacobian()
{
  return 0.;
}

Real TigErAntiDiffusionTerm::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}