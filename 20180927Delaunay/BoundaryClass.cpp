//*****************************************************************************
//Title		:Delaunay BoundaryClass
//Purpose	:Boundary for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/27
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************

#include "pch.h"
#include "BoundaryClass.h"


BoundaryClass::BoundaryClass(){}


BoundaryClass::~BoundaryClass(){}


//*****************************************************************************
//節点の境界における番号を返す
//*****************************************************************************
int BoundaryClass::order(int _nodenum) {
	for (int i = 0; i < nodelist.size(); i++) {
		if (nodelist[i] == _nodenum) {
			return i;
		}
	}
	return -1;
}