//*****************************************************************************
//Title		:src/cpp/Boundary.cpp
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#include "Boundary.h"


namespace DelaunayPAN{
	BoundaryClass::BoundaryClass(){}


	BoundaryClass::~BoundaryClass(){}


	//*****************************************************************************
	//�ߓ_�̋��E�ɂ�����ԍ���Ԃ�
	//*****************************************************************************
	int BoundaryClass::order(int _nodenum) {
		for (int i = 0; i < nodelist.size(); i++) {
			if (nodelist[i] == _nodenum) {
				return i;
			}
		}
		return -1;
	}
}