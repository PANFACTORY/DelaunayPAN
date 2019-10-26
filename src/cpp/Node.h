//*****************************************************************************
//Title		:src/cpp/Node.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>

namespace DelaunayPAN{
	class Node
	{
public:
		Node();
		~Node();

		double x, y;										//�ߓ_���W�_
		bool set;											//true�F�ݒu�ς�	false�F���ݒu
		bool type;											//true�F���E��̓_	false�F���E��̓_�łȂ�

		double distance(Node _node);					//�����̃m�[�h�Ƃ̋�����Ԃ�
		double vecpro(Node _node0, Node _node1);	//���g����_node0�ւ̃x�N�g����_node1�ւ̃x�N�g���̊O��
		double innpro(Node _node0, Node _node1);	//���g����_node0�ւ̃x�N�g����_node1�ւ̃x�N�g���̓���
	};
}