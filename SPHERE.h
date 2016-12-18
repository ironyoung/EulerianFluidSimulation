#ifndef SPHERE_H
#define SPHERE_H

#include "OBSTACLE.h"

class SPHERE : public OBSTACLE  
{
public:
	SPHERE(float x, float y, float z, float radius);
	virtual ~SPHERE();

  bool inside(float x, float y, float z);

private:
  float _center[3];
  float _radius;
};

#endif
