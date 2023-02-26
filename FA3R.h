#ifndef __FA3R_H__
#define __FA3R_H__

#pragma once

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/SVD"

using namespace std;
using namespace Eigen;

void FA3R_int(const vector<Vector3d>* P,
	          const vector<Vector3d>* Q,
	          Matrix3d * sigma,
              int num,
	          Matrix3d * rRes,
	          Vector3d * tRes);

void FA3R_double(const vector<Vector3d>* P,
	             const vector<Vector3d>* Q,
	             Matrix3d * sigma,
                 int num,
	             Matrix3d * rRes,
	             Vector3d * tRes);
	             
void eig3D_eig(const vector<Vector3d>* P,
	const vector<Vector3d>* Q,
	Matrix3d * sigma,
	Matrix3d * rRes,
	Vector3d * tRes);


#endif
