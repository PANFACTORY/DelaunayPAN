//*****************************************************************************
//Title		:src/cpp/Boundary.h
//Author	:Tanabe Yuta
//Date		:2019/10/26`
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


namespace DelaunayPAN{
	class Boundary
	{
public:
		Boundary();
		~Boundary();

		std::vector<int> nodelists;		//list of nodes on boundary
		bool type;						//type of boundary 

		int order(int _nodenum);		//get order id on boundary
	};


	Boundary::Boundary(){}


	Boundary::~Boundary(){}


	int Boundary::order(int _nodenum) {
		for (int i = 0; i < this->nodelists.size(); i++) {
			if (this->nodelists[i] == _nodenum) {
				return i;
			}
		}
		return -1;
	}
}