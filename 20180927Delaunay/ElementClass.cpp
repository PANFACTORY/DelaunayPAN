//*****************************************************************************
//Title		:Delaunay ElementClass
//Purpose	:Element for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/20
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************
#include "pch.h"
#include "ElementClass.h"
#include <float.h>
#include <cmath>


ElementClass::ElementClass(){
	side[0] = false;				//始め点は全て境界ではないものとする
	side[1] = false;
	side[2] = false;

	neighbor[0] = -1;				//隣接要素番号は始めは無し
	neighbor[1] = -1;
	neighbor[2] = -1;
}


ElementClass::~ElementClass(){}


//*****************************************************************************
//頂点を取り込む
//*****************************************************************************
void ElementClass::setnode(int _node0, int _node1, int _node2) {
	node[0] = _node0;
	node[1] = _node1;
	node[2] = _node2;
}


//*****************************************************************************
//隣接要素番号を取り込む
//*****************************************************************************
void ElementClass::setneighbor(int _neighbor0, int _neighbor1, int _neighbor2) {
	neighbor[0] = _neighbor0;
	neighbor[1] = _neighbor1;
	neighbor[2] = _neighbor2;
}


//*****************************************************************************
//境界辺を取り込む
//*****************************************************************************
void ElementClass::setside(bool _side0, bool _side1, bool _side2) {
	side[0] = _side0;
	side[1] = _side1;
	side[2] = _side2;
}


//*****************************************************************************
//引数のオブジェクトをコピーする
//*****************************************************************************
void ElementClass::copy(ElementClass _originalelement) {
	for (int i = 0; i < 3; i++) {
		node[i] = _originalelement.node[i];
		neighbor[i] = _originalelement.neighbor[i];
		side[i] = _originalelement.side[i];
	}
}


//*****************************************************************************
//三角形と節点の位置関係を判定
//	return -(i+1)	：i番目の辺の外側
//	return i+1		：i番目の辺上
//	return 0		：三角形内部
//*****************************************************************************

//**************************要修正*********************************************

int ElementClass::inouton(int _nodenum, vector<NodeClass> _node) {
	double vecpro0 = _node[node[0]].vecpro(_node[node[1]], _node[_nodenum]);
	double vecpro1 = _node[node[1]].vecpro(_node[node[2]], _node[_nodenum]);
	double vecpro2 = _node[node[2]].vecpro(_node[node[0]], _node[_nodenum]);

	//辺2外
	if (vecpro0 < 0.0 && vecpro2 > 0.0) {
		return -3;
	}
	//辺0外
	else if (vecpro1 < 0.0 && vecpro0 > 0.0) {
		return -1;
	}
	//辺1外
	else if (vecpro2 < 0.0 && vecpro1 > 0.0) {
		return -2;
	}
	//辺2上
	else if (fabs(vecpro0) <= DBL_EPSILON) {
		return 3;
	}
	//辺0上
	else if (fabs(vecpro1) <= DBL_EPSILON) {
		return 1;
	}
	//辺1上
	else if (fabs(vecpro2) <= DBL_EPSILON) {
		return 2;
	}
	return 0;
}


//*****************************************************************************
//指定された要素との位置関係を返す
//*****************************************************************************
int ElementClass::oppositenode(int _elementname) {
	for (int i = 0; i < 3; i++) {
		if (neighbor[i] == _elementname) {
			return i;
		}
	}
	return -1;
}


//*****************************************************************************
//頂角を計算
//*****************************************************************************
void ElementClass::getangle(vector<NodeClass> _node) {
	for (int i = 0; i < 3; i++) {
		angle[i] = _node[node[i]].vecpro(_node[node[(i + 1) % 3]], _node[node[(i + 2) % 3]]) / (_node[node[i]].distance(_node[node[(i + 1) % 3]]) * _node[node[i]].distance(_node[node[(i + 2) % 3]]));
	}
}