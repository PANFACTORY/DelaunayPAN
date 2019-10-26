//*****************************************************************************
//Title		:src/cpp/Node.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>


namespace DelaunayPAN{
	template<class T>
	class Node
	{
public:
		Node();
		~Node();

		T x, y;			//cordinate vallue of Node
		bool set;		//is set						
		bool type;		//is on boundary

		T distance(Node<T> _node);					//get distance between other node
		T vecpro(Node<T> _node0, Node<T> _node1);	//get innerproduct
		T innpro(Node<T> _node0, Node<T> _node1);	//get vectorproduct
	};


	template<class T>
	Node<T>::Node(){
		set = false;
		type = false;
	}


	template<class T>
	Node<T>::~Node(){}


	template<class T>
	T Node<T>::distance(Node<T> _node) {
		return sqrt(pow(x - _node.x, 2.0) + pow(y - _node.y, 2.0));
	}


	template<class T>
	T Node<T>::vecpro(Node<T> _node0, Node<T> _node1) {
		return (_node0.x - x)*(_node1.y - y) - (_node0.y - y)*(_node1.x - x);
	}


	template<class T>
	T Node<T>::innpro(Node<T> _node0, Node<T> _node1) {
		return (_node0.x - x)*(_node1.x - x) + (_node1.y - y)*(_node0.y - y);
	}
}