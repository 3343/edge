#ifndef IMPL_FXNS_H
#define IMPL_FXNS_H

#include "surf/Topo.h"

// TODO Incorporate global variables into a class
// NOTE These global variables are defined in "implicit_functions.h"
extern const double X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX, X0, Y0, F_ANGLE;
extern const std::string TOPO_FILE;
extern edge_cut::surf::Topo g_topo;


/**
 * Implicit function defining the region X_MIN <= x <= X_MAX, Y_MIN <= y <= Y_MAX,
 * Z_MIN <= z <= topography
 *
 * @param i_x x-coord of point to test
 * @param i_y y-coord of point to test
 * @param i_z z-coord of point to test
 * @return Negative value if point is inside region
 **/
double belowTopo( const double i_x,
                  const double i_y,
                  const double i_z );

/**
 * Implicit function defining the region X_MIN <= x <= X_MAX, faultPlane <= y <= Y_MAX,
 * Z_MIN <= z <= Z_MAX, where faultPlane is the surface (x-X_MIN)*sin(F_ANGLE) = y-Y0
 *
 * @param i_x x-coord of point to test
 * @param i_y y-coord of point to test
 * @param i_z z-coord of point to test
 * @return Negative value if point is inside region
 **/
double posHSpaceDisp( const double i_x,
                      const double i_y,
                      const double i_z  );

/**
 * Implicit function defining the region X_MIN <= x <= X_MAX, Y_MIN <= y <= faultPlane,
 * Z_MIN <= z <= Z_MAX, where faultPlane is the surface (x-X_MIN)*sin(F_ANGLE) = y-Y0
 *
 * @param i_x x-coord of point to test
 * @param i_y y-coord of point to test
 * @param i_z z-coord of point to test
 * @return Negative value if point is inside region
 **/
double negHSpaceDisp( const double i_x,
                      const double i_y,
                      const double i_z  );


/**
 * Template to convert standard function pointers into CGAL-compatible types
 **/
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
