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
	for (int i = 0; i < 3; i++) {
		double vecpro0 = _node[node[i % 3]].vecpro(_node[node[(i + 1) % 3]], _node[_nodenum]);
		double vecpro1 = _node[node[(i + 2) % 3]].vecpro(_node[node[i % 3]], _node[_nodenum]);
		
		if (fabs(vecpro0) <= DBL_EPSILON
			&& (_node[node[(i + 1) % 3]].x - _node[node[i % 3]].x)*(_node[_nodenum].x - _node[node[i % 3]].x) > 0.0
			&& fabs(_node[_nodenum].x - _node[node[i % 3]].x) < fabs(_node[node[(i + 1) % 3]].x - _node[node[i % 3]].x)) {	//辺上にあったとき
			return (i + 2) % 3 + 1;
		}

		else if (vecpro0 < 0.0 && vecpro1 > 0.0) {			//辺の外側にあったとき
			return -((i + 2) % 3 + 1);
		}
	}
	return 0;												//三角形内部にあったとき
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