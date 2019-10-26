//*****************************************************************************
//Title		:Delaunay NodeClass
//Purpose	:Node for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/20�`
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************

#include "NodeClass.h"
#include <cmath>


NodeClass::NodeClass(){
	set = false;			//���߂͓_�͐ݒu����Ă��Ȃ����̂Ƃ���
	type = false;			//���߂͋��E��̓_�ł͂Ȃ��Ƃ���
}


NodeClass::~NodeClass(){}


//*****************************************************************************
//�����Ɏw�肳�ꂽ�ߓ_�Ƃ̋�����Ԃ�
//*****************************************************************************
double NodeClass::distance(NodeClass _node) {
	return sqrt(pow(x - _node.x, 2.0) + pow(y - _node.y, 2.0));
}


//*****************************************************************************
//���g����_node0�ւ̃x�N�g����_node1�ւ̃x�N�g���̊O��
//*****************************************************************************
double NodeClass::vecpro(NodeClass _node0, NodeClass _node1) {
	return (_node0.x - x)*(_node1.y - y) - (_node0.y - y)*(_node1.x - x);
}


//*****************************************************************************
//���g����_node0�ւ̃x�N�g����_node1�ւ̃x�N�g���̊O��
//*****************************************************************************
double NodeClass::innpro(NodeClass _node0, NodeClass _node1) {
	return (_node0.x - x)*(_node1.x - x) + (_node1.y - y)*(_node0.y - y);
}
