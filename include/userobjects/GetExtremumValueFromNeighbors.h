#ifndef GETEXTREMUMVALUEFROMNEIGHBORS_H
#define GETEXTREMUMVALUEFROMNEIGHBORS_H

#include "NodalUserObject.h"

class GetExtremumValueFromNeighbors;
//class SymmTensor;

template<>
InputParameters validParams<GetExtremumValueFromNeighbors>();

class GetExtremumValueFromNeighbors : public NodalUserObject
{
public:
  GetExtremumValueFromNeighbors(const std::string & name, InputParameters parameters);

  ~GetExtremumValueFromNeighbors(); // the destructor closes the output file

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & u );
  virtual void finalize();

protected:

  // Auxiliary set of variables
  AuxiliarySystem & _aux;

  // Variable storing the nodal extremum value
  MooseVariable & _var;  
  MooseVariable & _var_out;

  // Boolean for maximum/minimum function
  bool _comp_max;

  // Real storing the extremum

  // The data structure used to find neighboring elements give a node ID
  std::vector< std::vector< const Elem * > > _nodes_to_elem_map;
};

#endif // GETEXTREMUMVALUEFROMNEIGHBORS_H




