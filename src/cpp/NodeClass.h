//*****************************************************************************
//Title		:Delaunay NodeClass
//Purpose	:Node for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/20`
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************
#pragma once
#include <vector>

class NodeClass
{
public:
	NodeClass();
	~NodeClass();

	double x, y;										//�ߓ_���W�_
	bool set;											//true�F�ݒu�ς�	false�F���ݒu
	bool type;											//true�F���E��̓_	false�F���E��̓_�łȂ�

	double distance(NodeClass _node);					//�����̃m�[�h�Ƃ̋�����Ԃ�
	double vecpro(NodeClass _node0, NodeClass _node1);	//���g����_node0�ւ̃x�N�g����_node1�ւ̃x�N�g���̊O��
	double innpro(NodeClass _node0, NodeClass _node1);	//���g����_node0�ւ̃x�N�g����_node1�ւ̃x�N�g���̓���
};

