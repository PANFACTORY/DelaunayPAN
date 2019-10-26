//*****************************************************************************
//Title		:src/cpp/Element.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#include "Element.h"
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>


namespace DelaunayPAN{
	#define SELFEPS0 DBL_EPSILON*1.0
	#define SELFEPS1 DBL_EPSILON*1.0


	Element::Element(){
		side[0] = false;				//�n�ߓ_�͑S�ċ��E�ł͂Ȃ����̂Ƃ���
		side[1] = false;
		side[2] = false;

		neighbor[0] = -1;				//�אڗv�f�ԍ��͎n�߂͖���
		neighbor[1] = -1;
		neighbor[2] = -1;

		active = true;
		check = false;
	}


	Element::~Element(){}


	//*****************************************************************************
	//���_����荞��
	//*****************************************************************************
	void Element::setnode(int _node0, int _node1, int _node2) {
		node[0] = _node0;
		node[1] = _node1;
		node[2] = _node2;
	}


	//*****************************************************************************
	//�אڗv�f�ԍ�����荞��
	//*****************************************************************************
	void Element::setneighbor(int _neighbor0, int _neighbor1, int _neighbor2) {
		neighbor[0] = _neighbor0;
		neighbor[1] = _neighbor1;
		neighbor[2] = _neighbor2;
	}


	//*****************************************************************************
	//���E�ӂ���荞��
	//*****************************************************************************
	void Element::setside(bool _side0, bool _side1, bool _side2) {
		side[0] = _side0;
		side[1] = _side1;
		side[2] = _side2;
	}


	//*****************************************************************************
	//�����̃I�u�W�F�N�g���R�s�[����
	//*****************************************************************************
	void Element::copy(Element _originalelement) {
		for (int i = 0; i < 3; i++) {
			node[i] = _originalelement.node[i];
			neighbor[i] = _originalelement.neighbor[i];
			side[i] = _originalelement.side[i];
			angle[i] = _originalelement.angle[i];
		}
	}


	//*****************************************************************************
	//�O�p�`�Ɛߓ_�̈ʒu�֌W�𔻒�
	//	return -(i+1)	�Fi�Ԗڂ̕ӂ̊O��
	//	return i+1		�Fi�Ԗڂ̕ӏ�
	//	return 0		�F�O�p�`����
	//*****************************************************************************
	int Element::inouton(int _nodenum, std::vector<Node> _node) {
		double vecpro0 = _node[node[0]].vecpro(_node[node[1]], _node[_nodenum]);
		double vecpro1 = _node[node[1]].vecpro(_node[node[2]], _node[_nodenum]);
		double vecpro2 = _node[node[2]].vecpro(_node[node[0]], _node[_nodenum]);

		double vecpro3 = _node[node[0]].vecpro(_node[node[2]], _node[_nodenum]);
		double vecpro4 = _node[node[1]].vecpro(_node[node[0]], _node[_nodenum]);
		double vecpro5 = _node[node[2]].vecpro(_node[node[1]], _node[_nodenum]);

		//��2�O
		if (vecpro2 > 0.0 && (vecpro0 < 0.0 || (fabs(vecpro0) <= SELFEPS0 && vecpro5 > 0.0))) {
			return -3;
		}
		//��0�O
		else if (vecpro0 > 0.0 && (vecpro1 < 0.0 || (fabs(vecpro1) <= SELFEPS0 && vecpro3 > 0.0))) {
			return -1;
		}
		//��1�O
		else if (vecpro1 > 0.0 && (vecpro2 < 0.0 || (fabs(vecpro2) <= SELFEPS0 && vecpro4 > 0.0))) {
			return -2;
		}
		//��2��
		else if (fabs(vecpro0) <= SELFEPS1) {
			return 3;
		}
		//��0��
		else if (fabs(vecpro1) <= SELFEPS1) {
			return 1;
		}
		//��1��
		else if (fabs(vecpro2) <= SELFEPS1) {
			return 2;
		}
		return 0;
	}


	//*****************************************************************************
	//�w�肳�ꂽ���_�����Ԗڂ��Ԃ�
	//*****************************************************************************
	int Element::nodeorder(int _nodenum) {
		for (int i = 0; i < 3; i++) {
			if (node[i] == _nodenum) {
				return i;
			}
		}
		return -1;
	}


	//*****************************************************************************
	//�w�肳�ꂽ�v�f�Ƃ̈ʒu�֌W��Ԃ�
	//*****************************************************************************
	int Element::oppositenode(int _elementname) {
		for (int i = 0; i < 3; i++) {
			if (neighbor[i] == _elementname) {
				return i;
			}
		}
		return -1;
	}


	//*****************************************************************************
	//���p���v�Z
	//*****************************************************************************
	void Element::getangle(std::vector<Node> _node) {
		for (int i = 0; i < 3; i++) {
			angle[i] = 180.0*acos(_node[node[i]].innpro(_node[node[(i + 1) % 3]], _node[node[(i + 2) % 3]]) / (_node[node[i]].distance(_node[node[(i + 1) % 3]]) * _node[node[i]].distance(_node[node[(i + 2) % 3]]))) / M_PI;
		}
	}


	//*****************************************************************************
	//�v�f�̖ʐς�Ԃ�
	//*****************************************************************************
	double Element::space(std::vector<Node> _node) {
		return 0.5*_node[node[0]].vecpro(_node[node[1]], _node[node[2]]);
	}
}