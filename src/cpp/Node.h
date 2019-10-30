//*****************************************************************************
//Title		:src/cpp/Node.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>


namespace DelaunayPAN{
	template<class T>
	class Node
	{
public:
		Node();
		~Node();
		Node(T _x, T _y);

		T x, y;					//cordinate vallue of Node
		bool isset;				//is set						
		bool isonboundary;		//is on boundary

		T distance(const Node<T>& _node);							//get distance between other node
		T vecpro(const Node<T>& _node0, const Node<T>& _node1);		//get innerproduct
		T innpro(const Node<T>& _node0, const Node<T>& _node1);		//get vectorproduct
	};


	template<class T>
	Node<T>::Node(){
		this->x = T();
		this->y = T();
		this->isset = false;
		this->isonboundary = false;
	}


	template<class T>
	Node<T>::~Node(){}


	template<class T>
	Node<T>::Node(T _x, T _y){
		this->x = _x;
		this->y = _y;
		this->isset = false;
		this->isonboundary = false;
	}


	template<class T>
	T Node<T>::distance(const Node<T>& _node) {
		return sqrt(pow(this->x - _node.x, 2.0) + pow(this->y - _node.y, 2.0));
	}


	template<class T>
	T Node<T>::vecpro(const Node<T>& _node0, const Node<T>& _node1) {
		return (_node0.x - this->x)*(_node1.y - this->y) - (_node0.y - this->y)*(_node1.x - this->x);
	}


	template<class T>
	T Node<T>::innpro(const Node<T>& _node0, const Node<T>& _node1) {
		return (_node0.x - this->x)*(_node1.x - this->x) + (_node1.y - this->y)*(_node0.y - this->y);
	}
}