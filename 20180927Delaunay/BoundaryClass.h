//*****************************************************************************
//Title		:Delaunay BoundaryClass
//Purpose	:Boundary for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/27～
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************

#include <vector>

using namespace std;

#pragma once
class BoundaryClass
{
public:
	BoundaryClass();
	~BoundaryClass();

	vector<int> nodelist;		//節点リスト
	bool type;					//境界の種類（true:外部　false:内部）

	int order(int _nodenum);	//節点の境界における番号を返す
};

