//*****************************************************************************
//Title		:src/cpp/Element.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <array>
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
		Element(int _node0, int _node1, int _node2);

		std::array<int, 3> nodes;								//id of nodes
		std::array<int, 3> neighbors;							//id of neighbor element
		std::array<T, 3> angles;								//angle of each corner
		std::array<bool, 3> sides;								//is edge on boundary
		bool active;											//is in boundary
		bool check;												//has already checked

		void setnode(int _node0, int _node1, int _node2);					//set node id
		void setneighbor(int _neighbor0, int _neighbor1, int _neighbor2);	//set neighbor id
		void setside(bool _side0, bool _side1, bool _side2);				//set edge is or not on boundary
		void copy(const Element<T>& _originalelement);						//copy element
		void getangle(std::vector<Node<T> >& _node);						//get angle

		int inouton(int _nodenum, std::vector<Node<T> >& _nodes);			//get location of node
		int oppositenode(int _elementname);									//get id of opposite node 

		T space(std::vector<Node<T> >& _node);								//get space of element
		int nodeorder(int _nodenum);										//get node id in element 
	};


	#define SELFEPS0 DBL_EPSILON*1.0
	#define SELFEPS1 DBL_EPSILON*1.0


	template<class T>
	Element<T>::Element(){
		this->nodes[0] = -1;
		this->nodes[1] = -1;
		this->nodes[2] = -1;

		this->sides[0] = false;
		this->sides[1] = false;
		this->sides[2] = false;

		this->neighbors[0] = -1;
		this->neighbors[1] = -1;
		this->neighbors[2] = -1;

		this->active = true;
		this->check = false;
	}


	template<class T>
	Element<T>::~Element(){}


	template<class T>
	Element<T>::Element(int _node0, int _node1, int _node2){
		this->nodes[0] = _node0;
		this->nodes[1] = _node1;
		this->nodes[2] = _node2;

		this->sides[0] = false;				
		this->sides[1] = false;
		this->sides[2] = false;

		this->neighbors[0] = -1;				
		this->neighbors[1] = -1;
		this->neighbors[2] = -1;

		this->active = true;
		this->check = false;
	}


	template<class T>
	void Element<T>::setnode(int _node0, int _node1, int _node2) {
		this->nodes[0] = _node0;
		this->nodes[1] = _node1;
		this->nodes[2] = _node2;
	}


	template<class T>
	void Element<T>::setneighbor(int _neighbor0, int _neighbor1, int _neighbor2) {
		this->neighbors[0] = _neighbor0;
		this->neighbors[1] = _neighbor1;
		this->neighbors[2] = _neighbor2;
	}


	template<class T>
	void Element<T>::setside(bool _side0, bool _side1, bool _side2) {
		this->sides[0] = _side0;
		this->sides[1] = _side1;
		this->sides[2] = _side2;
	}


	template<class T>
	void Element<T>::copy(const Element<T>& _originalelement) {
		this->nodes = _originalelement.nodes;
		this->neighbors = _originalelement.neighbors;
		this->sides = _originalelement.sides;
		this->angles = _originalelement.angles;
	}


	//*****************************************************************************
	//	return -(i+1)	�Fi�Ԗڂ̕ӂ̊O��
	//	return i+1		�Fi�Ԗڂ̕ӏ�
	//	return 0		�F�O�p�`����
	//*****************************************************************************
	template<class T>
	int Element<T>::inouton(int _nodenum, std::vector<Node<T> >& _nodes) {
		T vecpro0 = _nodes[this->nodes[0]].vecpro(_nodes[this->nodes[1]], _nodes[_nodenum]);
		T vecpro1 = _nodes[this->nodes[1]].vecpro(_nodes[this->nodes[2]], _nodes[_nodenum]);
		T vecpro2 = _nodes[this->nodes[2]].vecpro(_nodes[this->nodes[0]], _nodes[_nodenum]);

		T vecpro3 = _nodes[this->nodes[0]].vecpro(_nodes[this->nodes[2]], _nodes[_nodenum]);
		T vecpro4 = _nodes[this->nodes[1]].vecpro(_nodes[this->nodes[0]], _nodes[_nodenum]);
		T vecpro5 = _nodes[this->nodes[2]].vecpro(_nodes[this->nodes[1]], _nodes[_nodenum]);

		if (vecpro2 > T() && (vecpro0 < T() || (fabs(vecpro0) <= SELFEPS0 && vecpro5 > T()))) {
			return -3;
		}
		
		else if (vecpro0 > T() && (vecpro1 < T() || (fabs(vecpro1) <= SELFEPS0 && vecpro3 > T()))) {
			return -1;
		}
		
		else if (vecpro1 > T() && (vecpro2 < T() || (fabs(vecpro2) <= SELFEPS0 && vecpro4 > T()))) {
			return -2;
		}
		
		else if (fabs(vecpro0) <= SELFEPS1) {
			return 3;
		}
		
		else if (fabs(vecpro1) <= SELFEPS1) {
			return 1;
		}
		
		else if (fabs(vecpro2) <= SELFEPS1) {
			return 2;
		}
		return 0;
	}


	template<class T>
	int Element<T>::nodeorder(int _nodenum) {
		for (int i = 0; i < 3; i++) {
			if (this->nodes[i] == _nodenum) {
				return i;
			}
		}
		return -1;
	}


	template<class T>
	int Element<T>::oppositenode(int _elementname) {
		for (int i = 0; i < 3; i++) {
			if (this->neighbors[i] == _elementname) {
				return i;
			}
		}
		return -1;
	}


	template<class T>
	void Element<T>::getangle(std::vector<Node<T> >& _nodes) {
		for (int i = 0; i < 3; i++) {
			this->angles[i] = 180.0*acos(_nodes[this->nodes[i]].innpro(_nodes[this->nodes[(i + 1) % 3]], _nodes[this->nodes[(i + 2) % 3]]) / (_nodes[this->nodes[i]].distance(_nodes[this->nodes[(i + 1) % 3]]) * _nodes[this->nodes[i]].distance(_nodes[this->nodes[(i + 2) % 3]]))) / M_PI;
		}
	}


	template<class T>
	T Element<T>::space(std::vector<Node<T> >& _nodes) {
		return 0.5*_nodes[this->nodes[0]].vecpro(_nodes[this->nodes[1]], _nodes[this->nodes[2]]);
	}
}