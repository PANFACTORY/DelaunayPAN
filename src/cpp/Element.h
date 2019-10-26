//*****************************************************************************
//Title		:src/cpp/Element.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <float.h>
#define _USE_MATH_DEFINES
#include <cmath>


#include "Node.h"


namespace DelaunayPAN{
	template<class T>
	class Element{
public:
		Element();
		~Element();

		int node[3];											//id of nodes
		int neighbor[3];										//id of neighbor element
		T angle[3];												//angle of each corner
		bool side[3];											//is edge on boundary
		bool active;											//is in boundary
		bool check;												//has already checked

		void setnode(int _node0, int _node1, int _node2);					//set node id
		void setneighbor(int _neighbor0, int _neighbor1, int _neighbor2);	//set neighbor id
		void setside(bool _side0, bool _side1, bool _side2);				//set edge is or not on boundary
		void copy(Element<T> _originalelement);								//copy element
		void getangle(std::vector<Node<T> > _node);							//get angle

		int inouton(int _nodenum, std::vector<Node<T> > _node);				//get location of node
		int oppositenode(int _elementname);									//get id of opposite node 

		T space(std::vector<Node<T> > _node);								//get space of element
		int nodeorder(int _nodenum);										//get node id in element 
	};


	#define SELFEPS0 DBL_EPSILON*1.0
	#define SELFEPS1 DBL_EPSILON*1.0


	template<class T>
	Element<T>::Element(){
		side[0] = false;				//�n�ߓ_�͑S�ċ��E�ł͂Ȃ����̂Ƃ���
		side[1] = false;
		side[2] = false;

		neighbor[0] = -1;				//�אڗv�f�ԍ��͎n�߂͖���
		neighbor[1] = -1;
		neighbor[2] = -1;

		active = true;
		check = false;
	}


	template<class T>
	Element<T>::~Element(){}


	template<class T>
	void Element<T>::setnode(int _node0, int _node1, int _node2) {
		node[0] = _node0;
		node[1] = _node1;
		node[2] = _node2;
	}


	template<class T>
	void Element<T>::setneighbor(int _neighbor0, int _neighbor1, int _neighbor2) {
		neighbor[0] = _neighbor0;
		neighbor[1] = _neighbor1;
		neighbor[2] = _neighbor2;
	}


	template<class T>
	void Element<T>::setside(bool _side0, bool _side1, bool _side2) {
		side[0] = _side0;
		side[1] = _side1;
		side[2] = _side2;
	}


	template<class T>
	void Element<T>::copy(Element<T> _originalelement) {
		for (int i = 0; i < 3; i++) {
			node[i] = _originalelement.node[i];
			neighbor[i] = _originalelement.neighbor[i];
			side[i] = _originalelement.side[i];
			angle[i] = _originalelement.angle[i];
		}
	}


	//*****************************************************************************
	//	return -(i+1)	�Fi�Ԗڂ̕ӂ̊O��
	//	return i+1		�Fi�Ԗڂ̕ӏ�
	//	return 0		�F�O�p�`����
	//*****************************************************************************
	template<class T>
	int Element<T>::inouton(int _nodenum, std::vector<Node<T> > _node) {
		T vecpro0 = _node[node[0]].vecpro(_node[node[1]], _node[_nodenum]);
		T vecpro1 = _node[node[1]].vecpro(_node[node[2]], _node[_nodenum]);
		T vecpro2 = _node[node[2]].vecpro(_node[node[0]], _node[_nodenum]);

		T vecpro3 = _node[node[0]].vecpro(_node[node[2]], _node[_nodenum]);
		T vecpro4 = _node[node[1]].vecpro(_node[node[0]], _node[_nodenum]);
		T vecpro5 = _node[node[2]].vecpro(_node[node[1]], _node[_nodenum]);

		//��2�O
		if (vecpro2 > T() && (vecpro0 < T() || (fabs(vecpro0) <= SELFEPS0 && vecpro5 > T()))) {
			return -3;
		}
		//��0�O
		else if (vecpro0 > T() && (vecpro1 < T() || (fabs(vecpro1) <= SELFEPS0 && vecpro3 > T()))) {
			return -1;
		}
		//��1�O
		else if (vecpro1 > T() && (vecpro2 < T() || (fabs(vecpro2) <= SELFEPS0 && vecpro4 > T()))) {
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


	template<class T>
	int Element<T>::nodeorder(int _nodenum) {
		for (int i = 0; i < 3; i++) {
			if (node[i] == _nodenum) {
				return i;
			}
		}
		return -1;
	}


	template<class T>
	int Element<T>::oppositenode(int _elementname) {
		for (int i = 0; i < 3; i++) {
			if (neighbor[i] == _elementname) {
				return i;
			}
		}
		return -1;
	}


	template<class T>
	void Element<T>::getangle(std::vector<Node<T> > _node) {
		for (int i = 0; i < 3; i++) {
			angle[i] = 180.0*acos(_node[node[i]].innpro(_node[node[(i + 1) % 3]], _node[node[(i + 2) % 3]]) / (_node[node[i]].distance(_node[node[(i + 1) % 3]]) * _node[node[i]].distance(_node[node[(i + 2) % 3]]))) / M_PI;
		}
	}


	template<class T>
	T Element<T>::space(std::vector<Node<T> > _node) {
		return 0.5*_node[node[0]].vecpro(_node[node[1]], _node[node[2]]);
	}
}