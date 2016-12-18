#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <stdio.h>
#include <iostream>
#include <math.h>
#include "VEC3.h"
#include <assert.h>
// #include <VEC3.h>

namespace INTERPOLATE {

static inline float sameSign(float number, int sign)
{
	return (sign == 1) ? fabs(number) : (-fabs(number));
}

//////////////////////////////////////////////////////////////////////////////////////////
// cubic, bicubic, tricubic, n-cubic interpolation
// http://www.paulinternet.nl/?page=bicubic
//////////////////////////////////////////////////////////////////////////////////////////
static inline float monocubicInterpolate (float p[4], int x1, float x) {
	int	sign = 0;
	float dk0		= (p[2] - p[0])/2;
	float dk1		= (p[3] - p[1]) /2;
	float delta	= p[2] - p[1];
	if (delta == 0){
		dk0 = 0;
		dk1= 0;
	}
	else{
		sign	= (delta > 0) ? 1: -1;
		dk0		= sameSign(dk0, sign);
		dk1		= sameSign(dk1, sign);
	}
	const float a0	= p[1];
	const float a1	= dk0;
	const float a2	= 3 * delta - 2 * dk0 - dk1;
	const float a3	= dk0 + dk1 - delta;
	
	// x1 即为本位 
	return (a3*(x-x1)*(x-x1)*(x-x1) + a2*(x-x1)*(x-x1) + a1*(x-x1) + a0);
}

static inline float monobicubicInterpolate (float p[4][4], int x1, int y1, float x, float y) {
	float arr[4];
	arr[0] = monocubicInterpolate(p[0], y1, y);
	arr[1] = monocubicInterpolate(p[1], y1, y);
	arr[2] = monocubicInterpolate(p[2], y1, y);
	arr[3] = monocubicInterpolate(p[3], y1, y);
	return monocubicInterpolate(arr, x1, x);
}

static inline float monotricubicInterpolate (float p[4][4][4], int x1, int y1, int z1, float x, float y, float z) 
{
	float arr[4];
	arr[0] = monobicubicInterpolate(p[0], y1, z1, y, z);
	arr[1] = monobicubicInterpolate(p[1], y1, z1, y, z);
	arr[2] = monobicubicInterpolate(p[2], y1, z1, y, z);
	arr[3] = monobicubicInterpolate(p[3], y1, z1, y, z);
	return monocubicInterpolate(arr, x1, x);
}

static inline float monotriCubic3d(float* field, float x, float y, float z,  int xres, int yres, int zres)
{
	// clamp pos to grid boundaries
	if (x < 1.5) x = 1.5;
	if (x > xres - 2.5) x = xres - 2.5;
	if (y < 1.5) y = 1.5;
	if (y > yres - 2.5) y = yres - 2.5;
	if (z < 1.5) z = 1.5;
	if (z > zres - 2.5) z = zres - 2.5;

	// locate neighbors to interpolate
	const int x1 = (int)x;
	const int x0 = x1 - 1;
	const int x2 = x1 + 1;
	const int x3 = x2 + 1;
	const int xarray[4] = {x0, x1, x2, x3};

	const int y1 = (int)y;
	const int y0 = y1 - 1;
	const int y2 = y1 + 1;
	const int y3 = y2 + 1;
	const int yarray[4] = {y0, y1, y2, y3};

	const int z1 = (int)z;
	const int z0 = z1 - 1;
	const int z2 = z1 + 1;
	const int z3 = z2 + 1;
	const int zarray[4] = {z0, z1, z2, z3};

	const int slabSize	= xres*yres;
	int pos[4][4][4];
	float tmpfield[4][4][4];
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			for(int k = 0; k < 4; k++)
			{
				pos[i][j][k] = xarray[i] + yarray[j] * xres + zarray[k] * slabSize;
				tmpfield[i][j][k] = field[pos[i][j][k]];
			}

	return monotricubicInterpolate(tmpfield, x1, y1, z1, x, y, z);
}

static inline Vec3 monotriCubic3dVec(float* field1, float* field2, float* field3, 
		float x, float y, float z,  int xres, int yres, int zres)
{
	// clamp pos to grid boundaries
	if (x < 1.5) x = 1.5;
	if (x > xres - 2.5) x = xres - 2.5;
	if (y < 1.5) y = 1.5;
	if (y > yres - 2.5) y = yres - 2.5;
	if (z < 1.5) z = 1.5;
	if (z > zres - 2.5) z = zres - 2.5;

	// locate neighbors to interpolate
	const int x1 = (int)x;
	const int x0 = x1 - 1;
	const int x2 = x1 + 1;
	const int x3 = x2 + 1;
	const int xarray[4] = {x0, x1, x2, x3};

	const int y1 = (int)y;
	const int y0 = y1 - 1;
	const int y2 = y1 + 1;
	const int y3 = y2 + 1;
	const int yarray[4] = {y0, y1, y2, y3};

	const int z1 = (int)z;
	const int z0 = z1 - 1;
	const int z2 = z1 + 1;
	const int z3 = z2 + 1;
	const int zarray[4] = {z0, z1, z2, z3};

	const int slabSize	= xres*yres;
	int pos[4][4][4];
	float tmpfield1[4][4][4];
	float tmpfield2[4][4][4];
	float tmpfield3[4][4][4];
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			for(int k = 0; k < 4; k++)
			{
				pos[i][j][k] = xarray[i] + yarray[j] * xres + zarray[k] * slabSize;
				tmpfield1[i][j][k] = field1[pos[i][j][k]];
				tmpfield2[i][j][k] = field2[pos[i][j][k]];
				tmpfield3[i][j][k] = field3[pos[i][j][k]];
			}
	
	return Vec3(
			( monotricubicInterpolate(tmpfield1, x1, y1, z1, x, y, z) ) ,
			( monotricubicInterpolate(tmpfield2, x1, y1, z1, x, y, z) ) ,
			( monotricubicInterpolate(tmpfield3, x1, y1, z1, x, y, z) )
		);
}

//////////////////////////////////////////////////////////////////////
// linear interpolators
//////////////////////////////////////////////////////////////////////
static inline float lerp(float t, float a, float b) {
	return ( a + t * (b - a) );
}

// 2d linear interpolation
static inline float lerp(float* field, float x, float y, int res) {
	// clamp backtrace to grid boundaries
	if (x < 0.5f) x = 0.5f;
	if (x > res - 1.5f) x = res - 1.5f;
	if (y < 0.5f) y = 0.5f;
	if (y > res - 1.5f) y = res - 1.5f;

	const int x0 = (int)x;
	const int y0 = (int)y;
	x -= x0;
	y -= y0;
	float d00, d10, d01, d11;

	// lerp the velocities
	d00 = field[x0 + y0 * res];
	d10 = field[(x0 + 1) + y0 * res];
	d01 = field[x0 + (y0 + 1) * res];
	d11 = field[(x0 + 1) + (y0 + 1) * res];
	return lerp(y, lerp(x, d00, d10),
			lerp(x, d01, d11));
}

// 3d linear interpolation
static inline float lerp3d(float* field, float x, float y, float z,  int xres, int yres, int zres) {
	// clamp pos to grid boundaries
	if (x < 0.5) x = 0.5;
	if (x > xres - 1.5) x = xres - 1.5;
	if (y < 0.5) y = 0.5;
	if (y > yres - 1.5) y = yres - 1.5;
	if (z < 0.5) z = 0.5;
	if (z > zres - 1.5) z = zres - 1.5;

	// locate neighbors to interpolate
	const int x0 = (int)x;
	const int x1 = x0 + 1;
	const int y0 = (int)y;
	const int y1 = y0 + 1;
	const int z0 = (int)z;
	const int z1 = z0 + 1;

	// get interpolation weights
	const float s1 = x - (float)x0;
	const float s0 = 1.0f - s1;
	const float t1 = y - (float)y0;
	const float t0 = 1.0f - t1;
	const float u1 = z - (float)z0;
	const float u0 = 1.0f - u1;

	const int slabSize = xres*yres;
	const int i000 = x0 + y0 * xres + z0 * slabSize;
	const int i010 = x0 + y1 * xres + z0 * slabSize;
	const int i100 = x1 + y0 * xres + z0 * slabSize;
	const int i110 = x1 + y1 * xres + z0 * slabSize;
	const int i001 = x0 + y0 * xres + z1 * slabSize;
	const int i011 = x0 + y1 * xres + z1 * slabSize;
	const int i101 = x1 + y0 * xres + z1 * slabSize;
	const int i111 = x1 + y1 * xres + z1 * slabSize;

	// interpolate (indices could be computed once)
	return ( u0 * (s0 * (t0 * field[i000] +
		t1 * field[i010]) +
		s1 * (t0 * field[i100] +
		t1 * field[i110])) +
		u1 * (s0 * (t0 * field[i001] +
		t1 * field[i011]) +
		s1 * (t0 * field[i101] +
		t1 * field[i111])) );
}

// convert field entries of type T to floats, then interpolate
template <class T> 
static inline float lerp3dToFloat(T* field1,
		float x, float y, float z,  int xres, int yres, int zres) {
	// clamp pos to grid boundaries
	if (x < 0.5) x = 0.5;
	if (x > xres - 1.5) x = xres - 1.5;
	if (y < 0.5) y = 0.5;
	if (y > yres - 1.5) y = yres - 1.5;
	if (z < 0.5) z = 0.5;
	if (z > zres - 1.5) z = zres - 1.5;

	// locate neighbors to interpolate
	const int x0 = (int)x;
	const int x1 = x0 + 1;
	const int y0 = (int)y;
	const int y1 = y0 + 1;
	const int z0 = (int)z;
	const int z1 = z0 + 1;

	// get interpolation weights
	const float s1 = x - (float)x0;
	const float s0 = 1.0f - s1;
	const float t1 = y - (float)y0;
	const float t0 = 1.0f - t1;
	const float u1 = z - (float)z0;
	const float u0 = 1.0f - u1;

	const int slabSize = xres*yres;
	const int i000 = x0 + y0 * xres + z0 * slabSize;
	const int i010 = x0 + y1 * xres + z0 * slabSize;
	const int i100 = x1 + y0 * xres + z0 * slabSize;
	const int i110 = x1 + y1 * xres + z0 * slabSize;
	const int i001 = x0 + y0 * xres + z1 * slabSize;
	const int i011 = x0 + y1 * xres + z1 * slabSize;
	const int i101 = x1 + y0 * xres + z1 * slabSize;
	const int i111 = x1 + y1 * xres + z1 * slabSize;

	// interpolate (indices could be computed once)
	return (float)(
			( u0 * (s0 * (t0 * (float)field1[i000] +
				t1 * (float)field1[i010]) +
				s1 * (t0 * (float)field1[i100] +
				t1 * (float)field1[i110])) +
				u1 * (s0 * (t0 * (float)field1[i001] +
				t1 * (float)field1[i011]) +
				s1 * (t0 * (float)field1[i101] +
				t1 * (float)field1[i111])) ) );
}

static inline Vec3 lerp3dVec(float* field1, float* field2, float* field3, 
		float x, float y, float z,  int xres, int yres, int zres) {
	// clamp pos to grid boundaries
	if (x < 0.5) x = 0.5;
	if (x > xres - 1.5) x = xres - 1.5;
	if (y < 0.5) y = 0.5;
	if (y > yres - 1.5) y = yres - 1.5;
	if (z < 0.5) z = 0.5;
	if (z > zres - 1.5) z = zres - 1.5;

	// locate neighbors to interpolate
	const int x0 = (int)x;
	const int x1 = x0 + 1;
	const int y0 = (int)y;
	const int y1 = y0 + 1;
	const int z0 = (int)z;
	const int z1 = z0 + 1;

	// get interpolation weights
	const float s1 = x - (float)x0;
	const float s0 = 1.0f - s1;
	const float t1 = y - (float)y0;
	const float t0 = 1.0f - t1;
	const float u1 = z - (float)z0;
	const float u0 = 1.0f - u1;

	const int slabSize = xres*yres;
	const int i000 = x0 + y0 * xres + z0 * slabSize;
	const int i010 = x0 + y1 * xres + z0 * slabSize;
	const int i100 = x1 + y0 * xres + z0 * slabSize;
	const int i110 = x1 + y1 * xres + z0 * slabSize;
	const int i001 = x0 + y0 * xres + z1 * slabSize;
	const int i011 = x0 + y1 * xres + z1 * slabSize;
	const int i101 = x1 + y0 * xres + z1 * slabSize;
	const int i111 = x1 + y1 * xres + z1 * slabSize;

	// interpolate (indices could be computed once)
	return Vec3(
			( u0 * (s0 * (t0 * field1[i000] +
				t1 * field1[i010]) +
				s1 * (t0 * field1[i100] +
				t1 * field1[i110])) +
				u1 * (s0 * (t0 * field1[i001] +
				t1 * field1[i011]) +
				s1 * (t0 * field1[i101] +
				t1 * field1[i111])) ) , 
			( u0 * (s0 * (t0 * field2[i000] +
				t1 * field2[i010]) +
				s1 * (t0 * field2[i100] +
				t1 * field2[i110])) +
				u1 * (s0 * (t0 * field2[i001] +
				t1 * field2[i011]) +
				s1 * (t0 * field2[i101] +
				t1 * field2[i111])) ) , 
			( u0 * (s0 * (t0 * field3[i000] +
				t1 * field3[i010]) +
				s1 * (t0 * field3[i100] +
				t1 * field3[i110])) +
				u1 * (s0 * (t0 * field3[i001] +
				t1 * field3[i011]) +
				s1 * (t0 * field3[i101] +
				t1 * field3[i111])) ) 
			);
}

//////////////////////////////////////////////////////////////////////////////////////////
// hermite cubic interpolation
//////////////////////////////////////////////////////////////////////////////////////////
static inline Vec3 HermiteCubic3dVec(float* field1, float* field2, float* field3, 
		float x, float y, float z,  int xres, int yres, int zres) {
	// clamp pos to grid boundaries
	if (x < 1.5) x = 1.5;
	if (x > xres - 2.5) x = xres - 2.5;
	if (y < 1.5) y = 1.5;
	if (y > yres - 2.5) y = yres - 2.5;
	if (z < 1.5) z = 1.5;
	if (z > zres - 2.5) z = zres - 2.5;

	// locate neighbors to interpolate
	const int x0 = (int)x;
	const int x_1= x0 - 1;
	const int x1 = x0 + 1;
	const int x2 = x1 + 1;

	const int y0 = (int)y;
	const int y_1= y0 - 1;
	const int y1 = y0 + 1;
	const int y2 = y1 + 1;

	const int z0 = (int)z;
	const int z_1= z0 - 1;
	const int z1 = z0 + 1;
	const int z2 = z1 + 1;

	const int slabSize	= xres*yres;
	const int xpos0		= x0 + y0 * xres + z0 * slabSize;
	const int xpos_1	= x_1+ y0 * xres + z0 * slabSize;
	const int xpos1		= x1 + y0 * xres + z0 * slabSize;
	const int xpos2		= x2 + y0 * xres + z0 * slabSize;

	const int ypos0		= x0 + y0 * xres + z0 * slabSize;
	const int ypos_1	= x0 + y_1 * xres + z0 * slabSize;
	const int ypos1		= x0 + y1 * xres + z0 * slabSize;
	const int ypos2		= x0 + y2 * xres + z0 * slabSize;
	
	const int zpos0		= x0 + y0 * xres + z0 * slabSize;
	const int zpos_1	= x0 + y0 * xres + z_1* slabSize;
	const int zpos1		= x0 + y0 * xres + z1 * slabSize;
	const int zpos2		= x0 + y0 * xres + z2 * slabSize;

	// get interpolation coefficiences
	const float zdk		= (field3[zpos1] - field3[zpos_1])/2;
	const float zdk1	= (field3[zpos2] - field3[zpos0]) /2;
	const float zdelta	= field3[zpos1] - field3[zpos0];
	const float za0		= field3[zpos0];
	const float za1		= zdk;
	const float za2		= 3*zdelta - 2*zdk - zdk1;
	const float za3		= zdk + zdk1 - zdelta;

	const float ydk		= (field2[ypos1] - field2[ypos_1])/2;
	const float ydk1	= (field2[ypos2] - field2[ypos0])/2;
	const float ydelta	= field2[ypos1] - field2[ypos0];
	const float ya0		= field2[ypos0];
	const float ya1		= ydk;
	const float ya2		= 3*ydelta - 2*ydk - ydk1;
	const float ya3		= ydk + ydk1 - ydelta;

	const float xdk		= (field1[xpos1] - field1[xpos_1])/2;
	const float xdk1	= (field1[xpos2] - field1[xpos0])/2;
	const float xdelta	= field1[xpos1] - field1[xpos0];
	const float xa0		= field1[xpos0];
	const float xa1		= xdk;
	const float xa2		= 3*xdelta - 2*xdk - xdk1;
	const float xa3		= xdk + xdk1 - xdelta;

	// interpolate (indices could be computed once)
	return Vec3(
			( xa3*(x-x0)*(x-x0)*(x-x0) + xa2*(x-x0)*(x-x0) + xa1*(x-x0) + xa0 ) ,
			( ya3*(y-y0)*(y-y0)*(y-y0) + ya2*(y-y0)*(y-y0) + ya1*(y-y0) + ya0 ) ,
			( za3*(z-z0)*(z-z0)*(z-z0) + za2*(z-z0)*(z-z0) + za1*(z-z0) + za0 )
		);
}

//////////////////////////////////////////////////////////////////////////////////////////
// monotonic cubic interpolation
// < Visual Simulation of Smoke >
//////////////////////////////////////////////////////////////////////////////////////////
// 保证符号相等

//static inline Vec3 monoCubic3dVec(float* field1, float* field2, float* field3, 
//		float x, float y, float z,  int xres, int yres, int zres) {
//	// clamp pos to grid boundaries
//	if (x < 1.5) x = 1.5;
//	if (x > xres - 2.5) x = xres - 2.5;
//	if (y < 1.5) y = 1.5;
//	if (y > yres - 2.5) y = yres - 2.5;
//	if (z < 1.5) z = 1.5;
//	if (z > zres - 2.5) z = zres - 2.5;
//
//	// locate neighbors to interpolate
//	const int x0 = (int)x;
//	const int x_1= x0 - 1;
//	const int x1 = x0 + 1;
//	const int x2 = x1 + 1;
//
//	const int y0 = (int)y;
//	const int y_1= y0 - 1;
//	const int y1 = y0 + 1;
//	const int y2 = y1 + 1;
//
//	const int z0 = (int)z;
//	const int z_1= z0 - 1;
//	const int z1 = z0 + 1;
//	const int z2 = z1 + 1;
//
//	const int slabSize	= xres*yres;
//	const int xpos0	= x0 + y0 * xres + z0 * slabSize;
//	const int xpos_1= x_1+ y0 * xres + z0 * slabSize;
//	const int xpos1	= x1 + y0 * xres + z0 * slabSize;
//	const int xpos2	= x2 + y0 * xres + z0 * slabSize;
//
//	const int ypos0	= x0 + y0 * xres + z0 * slabSize;
//	const int ypos_1= x0 + y_1 * xres + z0 * slabSize;
//	const int ypos1	= x0 + y1 * xres + z0 * slabSize;
//	const int ypos2	= x0 + y2 * xres + z0 * slabSize;
//	
//	const int zpos0	= x0 + y0 * xres + z0 * slabSize;
//	const int zpos_1= x0 + y0 * xres + z_1* slabSize;
//	const int zpos1	= x0 + y0 * xres + z1 * slabSize;
//	const int zpos2	= x0 + y0 * xres + z2 * slabSize;
//
//	// get interpolation coefficiences
//	int	sign = 0;
//	float zdk		= (field3[zpos1] - field3[zpos_1])/2;
//	float zdk1		= (field3[zpos2] - field3[zpos0]) /2;
//	float zdelta	= field3[zpos1] - field3[zpos0];
//	if (zdelta == 0){
//		zdk = 0;
//		zdk1= 0;
//	}else{
//		sign	= (zdelta > 0) ? 1: -1;
//		zdk		= sameSign(zdk, sign);
//		zdk1	= sameSign(zdk1, sign);
//	}
//	const float za0	= field3[zpos0];
//	const float za1	= zdk;
//	const float za2	= 3*zdelta - 2*zdk - zdk1;
//	const float za3	= zdk + zdk1 - zdelta;
//
//	float ydk		= (field2[ypos1] - field2[ypos_1])/2;
//	float ydk1		= (field2[ypos2] - field2[ypos0])/2;
//	float ydelta	= field2[ypos1] - field2[ypos0];
//	if (ydelta == 0){
//		ydk = 0;
//		ydk1= 0;
//	}else{
//		sign	= (ydelta > 0) ? 1: -1;
//		ydk		= sameSign(ydk, sign);
//		ydk1	= sameSign(ydk1, sign);
//	}
//	const float ya0	= field2[ypos0];
//	const float ya1	= ydk;
//	const float ya2	= 3*ydelta - 2*ydk - ydk1;
//	const float ya3	= ydk + ydk1 - ydelta;
//
//	float xdk		= (field1[xpos1] - field1[xpos_1])/2;
//	float xdk1		= (field1[xpos2] - field1[xpos0])/2;
//	float xdelta	= field1[xpos1] - field1[xpos0];
//	if (xdelta == 0){
//		xdk = 0;
//		xdk1= 0;
//	}else{
//		sign	= (xdelta > 0) ? 1: -1;
//		xdk		= sameSign(xdk, sign);
//		xdk1	= sameSign(xdk1, sign);
//	}
//	const float xa0		= field1[xpos0];
//	const float xa1		= xdk;
//	const float xa2		= 3*xdelta - 2*xdk - xdk1;
//	const float xa3		= xdk + xdk1 - xdelta;
//
//	// interpolate (indices could be computed once)
//	return Vec3(
//			( xa3*(x-x0)*(x-x0)*(x-x0) + xa2*(x-x0)*(x-x0) + xa1*(x-x0) + xa0 ) ,
//			( ya3*(y-y0)*(y-y0)*(y-y0) + ya2*(y-y0)*(y-y0) + ya1*(y-y0) + ya0 ) ,
//			( za3*(z-z0)*(z-z0)*(z-z0) + za2*(z-z0)*(z-z0) + za1*(z-z0) + za0 )
//		);
//}

};
#endif
