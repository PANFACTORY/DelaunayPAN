//*****************************************************************************
//Title		:Delaunay ElementClass
//Purpose	:Element for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/20
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************
#pragma once
#include <vector>
#include "NodeClass.h"

using namespace std;

class ElementClass{
public:
	ElementClass();
	~ElementClass();

	int node[3];											//頂点の節点番号
	int neighbor[3];										//隣接要素番号
	bool side[3];											//true:辺が境界		false:辺が非境界

	void setnode(int _node0, int _node1, int _node2);		//頂点の取り込み
	void setneighbor(int _neighbor0, int _neighbor1, int _neighbor2);		//隣接要素番号の取り込み
	void setside(bool _side0, bool _side1, bool _side2);	//境界辺を取り込む
	void copy(ElementClass _originalelement);				//オブジェクトをコピーする

	int inouton(int _nodenum, vector<NodeClass> _node);		//要素と点の位置関係を返す
	int oppositenode(int _elementname);						//指定された要素との位置関係を返す
};

