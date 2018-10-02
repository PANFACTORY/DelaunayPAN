//*****************************************************************************
//Title		:Delaunay NodeClass
//Purpose	:Node for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/20
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************
#pragma once
#include <vector>

class NodeClass
{
public:
	NodeClass();
	~NodeClass();

	double x, y;										//節点座標点
	bool set;											//true：設置済み	false：未設置

	double distance(NodeClass _node);					//引数のノードとの距離を返す
	double vecpro(NodeClass _node0, NodeClass _node1);	//自身から_node0へのベクトルと_node1へのベクトルの外積
	double innpro(NodeClass _node0, NodeClass _node1);	//自身から_node0へのベクトルと_node1へのベクトルの内積
};

