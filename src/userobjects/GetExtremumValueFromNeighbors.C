#include "GetExtremumValueFromNeighbors.h"
//#include "SymmTensor.h"
//#include "FEProblem.h"
//#include <cmath>
//#include <algorithm>
//#include <set>

//libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/mesh_tools.h"

template<>
InputParameters validParams<GetExtremumValueFromNeighbors>()
{
  InputParameters params = validParams<NodalUserObject>();

  // Variables
  params.addRequiredParam<AuxVariableName>("variable", "Variable.");
  params.addRequiredParam<AuxVariableName>("variable_out", "Nodal variable storing the extremum value.");
  // Boolean for maximum/minimum function
  params.addParam<bool>("compute_maximum", true, "if true->maximum else false->minimum");

  return params;
}

GetExtremumValueFromNeighbors :: GetExtremumValueFromNeighbors(const std::string & name, InputParameters parameters) :
  NodalUserObject(name, parameters),
    // Auxiliary system
    _aux(_fe_problem.getAuxiliarySystem()),
    // Variables
    _var(_aux.getVariable(_tid, getParam<AuxVariableName>("variable"))),
    _var_out(_aux.getVariable(_tid, getParam<AuxVariableName>("variable_out"))),
    // Boolean for maximum/minimum function
    _comp_max(getParam<bool>("compute_maximum"))
{
}

GetExtremumValueFromNeighbors::~GetExtremumValueFromNeighbors()
{
}


void
GetExtremumValueFromNeighbors::initialize()
{
  // Initialyze to zero the vector storing 'variable'
  NumericVector<Number> & sln = _aux.solution();
  _aux.system().zero_variable(sln, _var_out.number());
  sln.close();

  // Build a new node to element map
  _nodes_to_elem_map.clear();
  MeshTools::build_nodes_to_elem_map(_mesh.getMesh(), _nodes_to_elem_map);
}

void
GetExtremumValueFromNeighbors::execute()
{
  // Get the neighbor nodes:
  std::vector< const Node * > neighbor_nodes;
  MeshTools::find_nodal_neighbors(_mesh.getMesh(), *_current_node, _nodes_to_elem_map, neighbor_nodes);

  // Loop over the neighbor nodes to determine the extremum value
  Real value = _is_implicit ? _var.getNodalValue(*_current_node) : _var.getNodalValueOld(*_current_node);
  for (unsigned int i_node=0; i_node<neighbor_nodes.size(); i_node++)
  {
    // Get the value at the node 'neighbor_node[i_node]'
    Number nghbr_nodal_val = _is_implicit ? _var.getNodalValue(*neighbor_nodes[i_node]) : _var.getNodalValueOld(*neighbor_nodes[i_node]);

    // If 'nodal_val' is larger/smaller than 'extremum', set 'extremum' equal to 'nodal_val'
    value = _comp_max ? std::max(nghbr_nodal_val,value) : std::min(nghbr_nodal_val,value);
  }

  // Store the computed extremum value 'value' in the variable called 'variable_out'
  NumericVector<Number> & sln = _aux.solution();  
  dof_id_type dof_out = _current_node->dof_number(_aux.number(), _var_out.number(), 0);
  sln.add(dof_out, value);
}


void
GetExtremumValueFromNeighbors::threadJoin(const UserObject & u )
{
}

void
GetExtremumValueFromNeighbors::finalize()
{
  _aux.solution().close();
}
