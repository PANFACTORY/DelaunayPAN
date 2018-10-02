//*****************************************************************************
//Title		:Delaunay NodeClass
//Purpose	:Node for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/20
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************
#include "pch.h"
#include "NodeClass.h"


NodeClass::NodeClass(){
	set = false;			//初めは点は設置されていないものとする
}


NodeClass::~NodeClass(){}


//*****************************************************************************
//引数に指定された節点との距離を返す
//*****************************************************************************
double NodeClass::distance(NodeClass _node) {
	return sqrt(pow(x - _node.x, 2.0) + pow(y - _node.y, 2.0));
}


//*****************************************************************************
//自身から_node0へのベクトルと_node1へのベクトルの外積
//*****************************************************************************
double NodeClass::vecpro(NodeClass _node0, NodeClass _node1) {
	return (_node0.x - x)*(_node1.y - y) - (_node0.y - y)*(_node1.x - x);
}


//*****************************************************************************
//自身から_node0へのベクトルと_node1へのベクトルの外積
//*****************************************************************************
double NodeClass::innpro(NodeClass _node0, NodeClass _node1) {
	return (_node0.x - x)*(_node1.x - x) + (_node1.y - y)*(_node0.y - y);
}
