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
This function computes the artificial dissipative terms. It only works in 1-D
 */
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
//    _kappa(getMaterialProperty<Real>("kappa")),
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
//  std::cout<<"&&&&&&&&&&&&&&"<<std::endl;
  // Get the nodal values and the shape functions
//  VariableValue var_nodal = _var.nodalValue();
  VariablePhiValue phi = _var.phi();
  VariablePhiGradient phi_grad = _var.gradPhi();
//  std::cout<<"phi values="<<phi[0][0]<<" and "<<phi[0][1]<<std::endl;
//  std::cout<<"grad phi values="<<phi_grad[0][0](0)<<" and "<<phi_grad[0][1](0)<<std::endl;
//  std::cout<<"c="<<_c<<std::endl;
//  std::cout<<"omega="<<_omega<<std::endl;
//  std::cout<<"JxW="<<_JxW[_qp]<<std::endl;  
//  std::cout<<"grad phi values="<<_grad_test[0][0](0)<<" and "<<_grad_test[0][1](0)<<std::endl;
//  std::cout<<"JxW="<<_JxW[_qp]<<std::endl;
//  std::cout<<var_nodal.size()<<std::endl;
//  std::cout<<phi.size()<<std::endl;
//  std::cout<<_phi.size()<<std::endl;
//  std::cout<<"nodal value="<<_u_nodal[0]<<" and "<<_u_nodal[1]<<std::endl;
//  std::cout<<"qp value="<<phi[0][0]*_u_nodal[0]+phi[1][0]*_u_nodal[1]<<" and "<<_u[0]<<std::endl;
//  std::cout<<"qp value="<<phi[0][1]*_u_nodal[0]+phi[1][1]*_u_nodal[1]<<" and "<<_u[1]<<std::endl;
//  std::cout<<"volume="<<_current_elem_volume<<std::endl;

  // Compute A_{i,j}, component (i,j) of the steady-state matrix
  std::vector<std::vector<Real> > Aij(phi.size(), std::vector<Real>(phi.size(), 0.));
//  Real Aij[phi.size()][phi.size()];
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
        Aij[ivar][jvar] += (_omega*phi_grad[jvar][qp](0)+_sigma[_qp]*phi[jvar][qp])*_c*_test[ivar][qp]*_coord[qp]*_JxW[qp];

  // Compute the bilinear form b_k(\phi_j, \phi_i)
  std::vector<std::vector<Real> > b_k(phi.size(), std::vector<Real>(phi.size(), 0.));
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      b_k[ivar][jvar] = jvar == ivar ? _current_elem_volume : -_current_elem_volume;

//  for (unsigned int jvar=0; jvar<phi.size(); jvar++)
//    for (unsigned int ivar=0; ivar<phi.size(); ivar++)
//    {
//      std::cout<<"Aij["<<ivar<<"]["<<jvar<<"]="<<Aij[ivar][jvar]<<std::endl;
//      std::cout<<"b_k["<<ivar<<"]["<<jvar<<"]="<<b_k[ivar][jvar]<<std::endl;
//    }

  // Compute the viscosity term
  Real sum_b_k(0.), max_Aij(0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      if (ivar != jvar)
      {
        max_Aij = std::max(max_Aij, std::max(0., Aij[ivar][jvar])/(-b_k[ivar][jvar]));
//        sum_b_k += b_k[ivar][jvar];
      }

//  std::cout<<"max_Aij="<<max_Aij<<std::endl;
//  std::cout<<"sum_b_k="<<sum_b_k<<std::endl;
  Real visc = max_Aij;//(-sum_b_k);
  
  if (visc<0)
    mooseError("The low-order viscosity coefficient computed in '"<<this->name()<<"' is locally negative.");
//  std::cout<<"visc="<<visc<<std::endl;

  // Compute the artificial diffusion term
//  Real diff_term_laplace = 0.5*_current_elem_volume*_grad_u[_qp](0)*_grad_test[_i][_qp](0);
//  std::cout<<"grad u="<<_grad_u[0](0)<<" and "<<_grad_u[1](0)<<std::endl;
//  std::cout<<"grad test="<<_grad_test[_i][_qp](0)<<std::endl;
  Real diff_term_laplace = _grad_u[_qp](0)*_grad_test[_i][_qp](0);
//  std::cout<<"diff_term_laplace="<<diff_term_laplace<<std::endl;

  Real diff_term(0.);
  for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    diff_term += _u_nodal[jvar]*b_k[_i][jvar];
  diff_term *= visc;//(_current_elem_volume*_current_elem_volume);
//  std::cout<<"diff_term="<<diff_term<<std::endl;

//  Real diff_term_la_total(0.), diff_term_total(0.);
//  for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
//    {
//      diff_term_la_total += _grad_u[qp](0)*_grad_test[_i][qp](0)*_JxW[qp]*_coord[qp];
//      diff_term_total += diff_term*_JxW[qp]*_coord[qp];
//    }

//  std::cout<<"diff_term="<<diff_term<<std::endl;
//  std::cout<<"diff_term_la total="<<diff_term_la_total*0.5*_current_elem_volume<<std::endl;

  // Return value
  return diff_term;
//  return diff_term_laplace*0.5*_current_elem_volume;
}

Real TigErArtificialVisc::computeQpJacobian()
{
  return 0.;
}

Real TigErArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}