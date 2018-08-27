#include "implicit_functions.h"

const double X_MIN =  -24950;
const double X_MAX =   24950;
const double Y_MIN =  -24950;
const double Y_MAX =   24950;
const double Z_MIN =  -10000;
const double Z_MAX =   80000;
const double X0    =  335000;
const double Y0    = 3775000;
const double F_ANGLE = 0;


double belowTopo( const double i_x,
                  const double i_y,
                  const double i_z )
{
  double l_xMinDisp   = X_MIN - i_x;
  double l_xMaxDisp   = i_x - X_MAX;
  double l_yMinDisp   = Y_MIN - i_y;
  double l_yMaxDisp   = i_y - Y_MAX;
  double l_zMinDisp   = Z_MIN - i_z;

  double l_max = std::max( {  l_xMinDisp,
                              l_xMaxDisp,
                              l_yMinDisp,
                              l_yMaxDisp,
                              l_zMinDisp  } );

  if ( l_max > 0 ) return l_max;
  else {
    double l_zTopoDisp  = g_topo.topoDisp( edge_cut::surf::TopoPoint( i_x, i_y, i_z ) );
    l_max = std::max( l_max, l_zTopoDisp );
  }

  return l_max;
}

double posHSpaceDisp( const double i_x,
                      const double i_y,
                      const double i_z  )
{
  double l_xMinDisp   = X_MIN - i_x;
  double l_xMaxDisp   = i_x - X_MAX;
  double l_yFaultDisp = Y0 + ( i_x - X0 ) * sin( F_ANGLE ) - i_y;
  double l_yMaxDisp   = i_y - Y_MAX;
  double l_zMinDisp   = Z_MIN - i_z;
  double l_zMaxDisp   = i_z - Z_MAX;

  double l_max = std::max( {  l_xMinDisp,
                              l_xMaxDisp,
                              l_yFaultDisp,
                              l_yMaxDisp,
                              l_zMinDisp,
                              l_zMaxDisp    } );
  return l_max;
}

double negHSpaceDisp( const double i_x,
                      const double i_y,
                      const double i_z  )
{
  double l_xMinDisp   = X_MIN - i_x;
  double l_xMaxDisp   = i_x - X_MAX;
  double l_yMinDisp   = Y_MIN - i_y;
  double l_yFaultDisp = i_y - (Y0 + ( i_x - X0 ) * sin( F_ANGLE ));
  double l_zMinDisp   = Z_MIN - i_z;
  double l_zMaxDisp   = i_z - Z_MAX;

  double l_max = std::max( {  l_xMinDisp,
                              l_xMaxDisp,
                              l_yMinDisp,
                              l_yFaultDisp,
                              l_zMinDisp,
                              l_zMaxDisp    } );
  return l_max;
}
