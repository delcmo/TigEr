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

#include "GetExtremumValueFromNeighbors.h"

template<>
InputParameters validParams<GetExtremumValueFromNeighbors>()
{
  InputParameters params = validParams<ElementUserObject>();

  // Variables
  params.addRequiredParam<NonlinearVariableName>("variable", "Variable.");
  params.addRequiredParam<AuxVariableName>("variable_out", "Nodal variable storing the extremum value.");
  // Boolean for maximum/minimum function
  params.addParam<bool>("compute_maximum", true, "if true->maximum else false->minimum");

  return params;
}

GetExtremumValueFromNeighbors::GetExtremumValueFromNeighbors(const std::string & name, InputParameters parameters) :
    ElementUserObject(name, parameters),
    // Set of nonlinear variables
    _nlsys(_fe_problem.getNonlinearSystem()),
    // Auxiliary system
    _aux(_fe_problem.getAuxiliarySystem()),
    // Variables
    _var(_nlsys.getVariable(_tid, parameters.get<NonlinearVariableName>("variable"))),
    _var_out(_aux.getVariable(_tid, getParam<AuxVariableName>("variable_out"))),
    // Boolean for maximum/minimum function
    _comp_max(getParam<bool>("compute_maximum"))
{
  mooseAssert(_mesh.dimension() != 1, "The function "<<_name<<" can only be used with a 1-D mesh");
}

void
GetExtremumValueFromNeighbors::initialize()
{
  // Initialyze to zero the vector storing 'variable'
  NumericVector<Number> & sln = _aux.solution();
  _aux.system().zero_variable(sln, _var_out.number());
  sln.close();
}

void
GetExtremumValueFromNeighbors::execute()
{
  mooseAssert(_current_elem->n_nodes() != 2, "The function "<<_name<<" can only be used with linear test function (two nodes per element).");

  // Get nodal values for nodes 'i' (0) and 'i+1' (1) belonging to '_current_elem'
  Number elem_nodal_val_i = _var.getNodalValue(*_current_elem->get_node(0));
  Number elem_nodal_val_ip1 = _var.getNodalValue(*_current_elem->get_node(1));
  Real extrem_value_elem = _comp_max ? std::max(elem_nodal_val_i, elem_nodal_val_ip1) : std::min(elem_nodal_val_i, elem_nodal_val_ip1);

  /// Compute extremum value for node 'i' (0) of '_current_element'
  // Determine neighbor element for node 'i' (left)
  const Elem * nghb_elem_to_node_i = _current_elem->neighbor(0) == NULL ? _current_elem : _current_elem->neighbor(0);

  // Get the nodal values 'i-1' and 'i' belonging to 'nghb_elem_to_node_i'
  Number nghb_nodal_val_im1 = _var.getNodalValue(*nghb_elem_to_node_i->get_node(0));
  Number nghb_nodal_val_i = _var.getNodalValue(*nghb_elem_to_node_i->get_node(1));

  // Determine extremum value for node 'i' of '_current_elem'
  Real extrem_value_nghb = _comp_max ? std::max(nghb_nodal_val_i, nghb_nodal_val_im1) : std::min(nghb_nodal_val_i, nghb_nodal_val_im1);
  Real extrem_value = _comp_max ? std::max(extrem_value_nghb, extrem_value_elem) : std::min(extrem_value_nghb, extrem_value_elem);

  // Store the computed extremum value 'extrem_value' in the variable called 'variable_out'
  NumericVector<Number> & sln = _aux.solution();
  dof_id_type dof_out_i = _current_elem->get_node(0)->dof_number(_aux.number(), _var_out.number(), 0);
  sln.set(dof_out_i, extrem_value);

  /// Compute extremum value for node 'i+1' (1) of '_current_element' only if it is last element of the mesh:
  if (_current_elem->neighbor(1) == NULL)
  {
    // Store the computed extremum value 'extrem_value' in the variable called 'variable_out'
    dof_id_type dof_out_ip1 = _current_elem->get_node(1)->dof_number(_aux.number(), _var_out.number(), 0);
    sln.set(dof_out_ip1, extrem_value_elem);
  }
}

void
GetExtremumValueFromNeighbors::finalize()
{
  _aux.solution().close();  
}