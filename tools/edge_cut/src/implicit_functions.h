#ifndef IMPL_FXNS_H
#define IMPL_FXNS_H

#include "surf/Topo.h"

extern const double X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX, X0, Y0, F_ANGLE;
extern const std::string TOPO_FILE;
extern edge_cut::surf::Topo g_topo;

// double topoDisp(  const edge_cut::surf::Topo &  i_topo,
//                   const double                  i_x,
//                   const double                  i_y,
//                   const double                  i_z     );

double belowTopo( const double i_x,
                  const double i_y,
                  const double i_z );

double posHSpaceDisp( const double i_x,
                      const double i_y,
                      const double i_z  );

double negHSpaceDisp( const double i_x,
                      const double i_y,
                      const double i_z  );

template <typename FT, typename P>
class FT_to_point_function_wrapper : public std::unary_function<P, FT>
{
  typedef FT (*Implicit_function)(FT, FT, FT);
  Implicit_function function;
public:
  typedef P Point;
  FT_to_point_function_wrapper(Implicit_function f) : function(f) {}
  FT operator()(Point p) const { return function(p.x(), p.y(), p.z()); }
};

#endif
