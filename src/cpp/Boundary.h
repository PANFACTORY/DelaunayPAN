//*****************************************************************************
//Title		:src/cpp/Boundary.h
//Author	:Tanabe Yuta
//Date		:2019/10/26`
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


namespace DelaunayPAN{
	class BoundaryClass
	{
public:
		BoundaryClass();
		~BoundaryClass();

		std::vector<int> nodelist;		//list of nodes on boundary
		bool type;						//type of boundary 

		int order(int _nodenum);		//get order id on boundary
	};
}