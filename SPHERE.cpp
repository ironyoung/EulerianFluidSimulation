#include "SPHERE.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SPHERE::SPHERE(float x, float y, float z, float radius) :
  _radius(radius)
{
  _center[0] = x;
  _center[1] = y;
  _center[2] = z;
}

SPHERE::~SPHERE()
{

}

bool SPHERE::inside(float x, float y, float z)
{
  float translate[] = {x - _center[0], y - _center[1], z - _center[2]};
  // float magnitude = translate[0] * translate[0] + 
  //                   translate[1] * translate[1] + 
  //                   translate[2] * translate[2];

    float magnitude = translate[0] * translate[0] + 
                    translate[1] * translate[1];

  return (magnitude < _radius * _radius) ? true : false;
}
