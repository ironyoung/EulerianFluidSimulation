#include "FLUID_3D.h"
#define SOLVER_ACCURACY 1e-06

//////////////////////////////////////////////////////////////////////
// solve the poisson equation with CG
//////////////////////////////////////////////////////////////////////
void FLUID_3D::solvePressure(float* field, float* b, unsigned char* skip)
{
  int x, y, z, index;

  // i = 0
  int i = 0;

  // r = b - Ax
  index = _slabSize + _xRes + 1;
  for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
      {
        // if the cell is a variable
        float Acenter = 0.0f;
        if (!skip[index])
        {
          // set the matrix to the Poisson stencil in order
          if (!skip[index + 1]) Acenter += 1.;
          if (!skip[index - 1]) Acenter += 1.;
          if (!skip[index + _xRes]) Acenter += 1.;
          if (!skip[index - _xRes]) Acenter += 1.;
          if (!skip[index + _slabSize]) Acenter += 1.;
          if (!skip[index - _slabSize]) Acenter += 1.;
        }
        
        _residual[index] = b[index] - (Acenter * field[index] +  
          field[index - 1] * (skip[index - 1] ? 0.0 : -1.0f)+ 
          field[index + 1] * (skip[index + 1] ? 0.0 : -1.0f)+
          field[index - _xRes] * (skip[index - _xRes] ? 0.0 : -1.0f)+ 
          field[index + _xRes] * (skip[index + _xRes] ? 0.0 : -1.0f)+
          field[index - _slabSize] * (skip[index - _slabSize] ? 0.0 : -1.0f)+ 
          field[index + _slabSize] * (skip[index + _slabSize] ? 0.0 : -1.0f) );
        _residual[index] = (skip[index]) ? 0.0f : _residual[index];
      }

  // d = r
  index = _slabSize + _xRes + 1;
  for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
        _direction[index] = _residual[index];

  // deltaNew = transpose(r) * r
  float deltaNew = 0.0f;
  index = _slabSize + _xRes + 1;
  for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
        deltaNew += _residual[index] * _residual[index];

  // delta0 = deltaNew
  float delta0 = deltaNew;

  // While deltaNew > (eps^2) * delta0
  const float eps  = SOLVER_ACCURACY;
  float maxR = 2.0f * eps;
  while ((i < _iterations) && (maxR > eps))
  {
    // q = Ad
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
        {
          // if the cell is a variable
          float Acenter = 0.0f;
          if (!skip[index])
          {
            // set the matrix to the Poisson stencil in order
            if (!skip[index + 1]) Acenter += 1.;
            if (!skip[index - 1]) Acenter += 1.;
            if (!skip[index + _xRes]) Acenter += 1.;
            if (!skip[index - _xRes]) Acenter += 1.;
            if (!skip[index + _slabSize]) Acenter += 1.;
            if (!skip[index - _slabSize]) Acenter += 1.;
          }
          
          _q[index] = Acenter * _direction[index] +  
            _direction[index - 1] * (skip[index - 1] ? 0.0 : -1.0f) + 
            _direction[index + 1] * (skip[index + 1] ? 0.0 : -1.0f) +
            _direction[index - _xRes] * (skip[index - _xRes] ? 0.0 : -1.0f) + 
            _direction[index + _xRes] * (skip[index + _xRes] ? 0.0 : -1.0f)+
            _direction[index - _slabSize] * (skip[index - _slabSize] ? 0.0 : -1.0f) + 
            _direction[index + _slabSize] * (skip[index + _slabSize] ? 0.0 : -1.0f);
          _q[index] = (skip[index]) ? 0.0f : _q[index];
        }

    // alpha = deltaNew / (transpose(d) * q)
    float alpha = 0.0f;
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          alpha += _direction[index] * _q[index];
    if (fabs(alpha) > 0.0f)
      alpha = deltaNew / alpha;

    // x = x + alpha * d
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          field[index] += alpha * _direction[index];

    // r = r - alpha * q
    maxR = 0.0f;
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
        {
          _residual[index] -= alpha * _q[index];
          maxR = (_residual[index] > maxR) ? _residual[index] : maxR;
        }

    // deltaOld = deltaNew
    float deltaOld = deltaNew;

    // deltaNew = transpose(r) * r
    deltaNew = 0.0f;
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          deltaNew += _residual[index] * _residual[index];

    // beta = deltaNew / deltaOld
    float beta = deltaNew / deltaOld;

    // d = r + beta * d
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          _direction[index] = _residual[index] + beta * _direction[index];

    // i = i + 1
    i++;
  }
  cout << i << " iterations converged to " << maxR << endl;
  if(i >= _iterations)
      system("pause");
}

//////////////////////////////////////////////////////////////////////
// solve the heat equation with CG
//////////////////////////////////////////////////////////////////////
void FLUID_3D::solveHeat(float* field, float* b, unsigned char* skip)
{
  int x, y, z, index;
  const float heatConst = _dt * _heatDiffusion / (_dx * _dx);

  // i = 0
  int i = 0;

  // r = b - Ax
  index = _slabSize + _xRes + 1;
  for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
      {
        // if the cell is a variable
        float Acenter = 1.0f;
        if (!skip[index])
        {
          // set the matrix to the Poisson stencil in order
          if (!skip[index + 1]) Acenter += heatConst;
          if (!skip[index - 1]) Acenter += heatConst;
          if (!skip[index + _xRes]) Acenter += heatConst;
          if (!skip[index - _xRes]) Acenter += heatConst;
          if (!skip[index + _slabSize]) Acenter += heatConst;
          if (!skip[index - _slabSize]) Acenter += heatConst;
        }
        
        _residual[index] = b[index] - (Acenter * field[index] + 
          field[index - 1] * (skip[index - 1] ? 0.0 : -heatConst) + 
          field[index + 1] * (skip[index + 1] ? 0.0 : -heatConst) +
          field[index - _xRes] * (skip[index - _xRes] ? 0.0 : -heatConst) + 
          field[index + _xRes] * (skip[index + _xRes] ? 0.0 : -heatConst) +
          field[index - _slabSize] * (skip[index - _slabSize] ? 0.0 : -heatConst) + 
          field[index + _slabSize] * (skip[index + _slabSize] ? 0.0 : -heatConst));
        _residual[index] = (skip[index]) ? 0.0f : _residual[index];
      }

  // d = r
  index = _slabSize + _xRes + 1;
  for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
        _direction[index] = _residual[index];

  // deltaNew = transpose(r) * r
  float deltaNew = 0.0f;
  index = _slabSize + _xRes + 1;
  for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
        deltaNew += _residual[index] * _residual[index];

  // delta0 = deltaNew
  float delta0 = deltaNew;

  // While deltaNew > (eps^2) * delta0
  const float eps  = SOLVER_ACCURACY;
  float maxR = 2.0f * eps;
  while ((i < _iterations) && (maxR > eps))
  {
    // q = Ad
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
        {
          // if the cell is a variable
          float Acenter = 1.0f;
          if (!skip[index])
          {
            // set the matrix to the Poisson stencil in order
            if (!skip[index + 1]) Acenter += heatConst;
            if (!skip[index - 1]) Acenter += heatConst;
            if (!skip[index + _xRes]) Acenter += heatConst;
            if (!skip[index - _xRes]) Acenter += heatConst;
            if (!skip[index + _slabSize]) Acenter += heatConst;
            if (!skip[index - _slabSize]) Acenter += heatConst;
          }

          _q[index] = (Acenter * _direction[index] + 
            _direction[index - 1] * (skip[index - 1] ? 0.0 : -heatConst) + 
            _direction[index + 1] * (skip[index + 1] ? 0.0 : -heatConst) +
            _direction[index - _xRes] * (skip[index - _xRes] ? 0.0 : -heatConst) + 
            _direction[index + _xRes] * (skip[index + _xRes] ? 0.0 : -heatConst) +
            _direction[index - _slabSize] * (skip[index - _slabSize] ? 0.0 : -heatConst) + 
            _direction[index + _slabSize] * (skip[index + _slabSize] ? 0.0 : -heatConst));
         
          _q[index] = (skip[index]) ? 0.0f : _q[index];
        }

    // alpha = deltaNew / (transpose(d) * q)
    float alpha = 0.0f;
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          alpha += _direction[index] * _q[index];
    if (fabs(alpha) > 0.0f)
      alpha = deltaNew / alpha;

    // x = x + alpha * d
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          field[index] += alpha * _direction[index];

    // r = r - alpha * q
    maxR = 0.0f;
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
        {
          _residual[index] -= alpha * _q[index];
          maxR = (_residual[index] > maxR) ? _residual[index] : maxR;
        }

    // deltaOld = deltaNew
    float deltaOld = deltaNew;

    // deltaNew = transpose(r) * r
    deltaNew = 0.0f;
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          deltaNew += _residual[index] * _residual[index];

    // beta = deltaNew / deltaOld
    float beta = deltaNew / deltaOld;

    // d = r + beta * d
    index = _slabSize + _xRes + 1;
    for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
      for (y = 1; y < _yRes - 1; y++, index += 2)
        for (x = 1; x < _xRes - 1; x++, index++)
          _direction[index] = _residual[index] + beta * _direction[index];

    // i = i + 1
    i++;
  }
  cout << i << " iterations converged to " << maxR << endl;
}
