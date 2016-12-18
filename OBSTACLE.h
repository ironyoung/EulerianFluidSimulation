#ifndef OBSTACLE_H
#define OBSTACLE_H

enum OBSTACLE_FLAGS {
	EMPTY = 0, 
	MARCHED = 2, 
	RETIRED = 4 
};  

class OBSTACLE  
{
public:
	OBSTACLE() {};
	virtual ~OBSTACLE() {};

  virtual bool inside(float x, float y, float z) = 0;
};

#endif
