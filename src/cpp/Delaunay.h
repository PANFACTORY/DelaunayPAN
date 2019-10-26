//*****************************************************************************
//Title		:src/cpp/Delaunay.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>


#include "Node.h"
#include "Element.h"
#include "Boundary.h"


namespace DelaunayPAN{
	template<class T>
	class Delaunay
	{
public:
		Delaunay();
		~Delaunay();

		void delaunaymain(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, std::vector<Boundary> &_boundaries, T _maxsize, int _laplaciannum);

		void getsupertriangle(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements);
		void deletesupertriangle(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element);

		void getboundary(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, Boundary _boundary);
		void deactivate(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, Boundary _boundary);
		void deleteelement(std::vector<Element<T> > &_elements);
		void sortelement(std::vector<Element<T> > &_elements);
		void getinternalelement(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, T _maxside);
		void swapping(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, std::vector<int> &_stack, int _nodenump1, int _nodenumm1);

		void getelementin(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1);
		void getelementon(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1);

		void laplacian(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, int _maxnum);
	};


	#define ADDITIONALNODENUM0	100000			//�ӂ̒����ɂ��ו�����


	template<class T>
	Delaunay<T>::Delaunay() {}


	template<class T>
	Delaunay<T>::~Delaunay() {}


	template<class T>
	void Delaunay<T>::delaunaymain(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, std::vector<Boundary> &_boundaries, T _maxsize, int _laplaciannum) {
		//----------Generate SuperTriangle----------
		getsupertriangle(_nodes, _elements);

		//----------Generate Boundary----------
		for (auto& boundary : _boundaries) {
			getboundary(_nodes, _elements, boundary);
		}
		for (auto& boundary : _boundaries) {
			deactivate(_nodes, _elements, boundary);
		}

		//----------Delete needless Elements----------
		deletesupertriangle(_nodes, _elements);
		deleteelement(_elements);

		//----------Sort Elements----------
		sortelement(_elements);

		if (_maxsize > 0) {
			//----------subdivide----------
			getinternalelement(_nodes, _elements, _maxsize);

			//----------Laplacian smoothing----------
			laplacian(_nodes, _elements, _laplaciannum);
		}
	}


	template<class T>
	void Delaunay<T>::getelementin(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1) {
		std::vector<int> stack;
		std::vector<Element<T> > tmpelement(3);

		tmpelement[0].sides[0] = _elements[_nowtri].sides[0];
		tmpelement[1].sides[0] = _elements[_nowtri].sides[1];
		tmpelement[2].sides[0] = _elements[_nowtri].sides[2];

		//���E�Ӕ���
		if (_elements[_nowtri].nodes[0] == _nodenumm1 || _elements[_nowtri].nodes[0] == _nodenump1) {
			tmpelement[1].sides[1] = true;
			tmpelement[2].sides[2] = true;
		}
		if (_elements[_nowtri].nodes[1] == _nodenumm1 || _elements[_nowtri].nodes[1] == _nodenump1) {
			tmpelement[2].sides[1] = true;
			tmpelement[0].sides[2] = true;
		}
		if (_elements[_nowtri].nodes[2] == _nodenumm1 || _elements[_nowtri].nodes[2] == _nodenump1) {
			tmpelement[0].sides[1] = true;
			tmpelement[1].sides[2] = true;
		}

		tmpelement[0].setnode(_nodenum, _elements[_nowtri].nodes[1], _elements[_nowtri].nodes[2]);
		tmpelement[0].setneighbor(_elements[_nowtri].neighbors[0], _elements.size(), _elements.size() + 1);

		tmpelement[1].setnode(_nodenum, _elements[_nowtri].nodes[2], _elements[_nowtri].nodes[0]);
		tmpelement[1].setneighbor(_elements[_nowtri].neighbors[1], _elements.size() + 1, _nowtri);

		tmpelement[2].setnode(_nodenum, _elements[_nowtri].nodes[0], _elements[_nowtri].nodes[1]);
		tmpelement[2].setneighbor(_elements[_nowtri].neighbors[2], _nowtri, _elements.size());

		tmpelement[0].getangle(_nodes);
		tmpelement[1].getangle(_nodes);
		tmpelement[2].getangle(_nodes);

		for (int k = 0; k < 2; k++) {
			int neighbor = _elements[_nowtri].neighbors[1 + k];
			if (neighbor >= 0) {
				_elements[neighbor].neighbors[_elements[neighbor].oppositenode(_nowtri)] = _elements.size() + k;
			}
		}

		stack.push_back(_nowtri);
		stack.push_back(_elements.size());
		stack.push_back(_elements.size() + 1);

		_elements[_nowtri].copy(tmpelement[0]);
		_elements.push_back(tmpelement[1]);
		_elements.push_back(tmpelement[2]);

		//.....�X���b�s���O.....
		swapping(_nodes, _elements, stack, _nodenump1, _nodenumm1);
	}


	template<class T>
	void Delaunay<T>::getelementon(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element, int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1) {
		std::vector<int> stack;
		int nownode = _pos;
		int neitri = _element[_nowtri].neighbors[nownode];
		//�ӂ�����ŗאڗv�f�����݂���Ƃ�
		if (neitri != -1 && _element[neitri].active == true) {
			int neinode = _element[neitri].oppositenode(_nowtri);
			std::vector<Element<T> > tmptri(4);			//0:nowtri	1:neitri

			tmptri[0].setnode(_nodenum, _element[_nowtri].nodes[nownode], _element[_nowtri].nodes[(nownode + 1) % 3]);
			tmptri[0].setneighbor(_element[_nowtri].neighbors[(nownode + 2) % 3], neitri, _element.size());
			tmptri[0].setside(_element[_nowtri].sides[(nownode + 2) % 3], false, false);

			tmptri[1].setnode(_nodenum, _element[_nowtri].nodes[(nownode + 1) % 3], _element[neitri].nodes[neinode]);
			tmptri[1].setneighbor(_element[neitri].neighbors[(neinode + 1) % 3], _element.size() + 1, _nowtri);
			tmptri[1].setside(_element[neitri].sides[(neinode + 1) % 3], false, false);

			tmptri[2].setnode(_nodenum, _element[_nowtri].nodes[(nownode + 2) % 3], _element[_nowtri].nodes[nownode]);
			tmptri[2].setneighbor(_element[_nowtri].neighbors[(nownode + 1) % 3], _nowtri, _element.size() + 1);
			tmptri[2].sides[0] = _element[_nowtri].sides[(nownode + 1) % 3];

			tmptri[3].setnode(_nodenum, _element[neitri].nodes[neinode], _element[_nowtri].nodes[(nownode + 2) % 3]);
			tmptri[3].setneighbor(_element[neitri].neighbors[(neinode + 2) % 3], _element.size(), neitri);
			tmptri[3].sides[0] = _element[neitri].sides[(neinode + 2) % 3];

			int nei1 = _element[_nowtri].neighbors[(nownode + 1) % 3];
			if (nei1 != -1) {
				_element[nei1].neighbors[_element[nei1].oppositenode(_nowtri)] = _element.size();
			}

			int nei2 = _element[neitri].neighbors[(neinode + 2) % 3];
			if (nei2 != -1) {
				_element[nei2].neighbors[_element[nei2].oppositenode(neitri)] = _element.size() + 1;
			}

			if (tmptri[0].nodes[1] == _nodenumm1 || tmptri[0].nodes[1] == _nodenump1) {
				tmptri[2].sides[1] = true;
				tmptri[0].sides[2] = true;
			}
			if (tmptri[1].nodes[1] == _nodenumm1 || tmptri[1].nodes[1] == _nodenump1) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
			}
			if (tmptri[2].nodes[1] == _nodenumm1 || tmptri[2].nodes[1] == _nodenump1) {
				tmptri[2].sides[2] = true;
				tmptri[3].sides[1] = true;
			}
			if (tmptri[3].nodes[1] == _nodenumm1 || tmptri[3].nodes[1] == _nodenump1) {
				tmptri[3].sides[2] = true;
				tmptri[1].sides[1] = true;
			}

			if (_element[_nowtri].sides[nownode] == true) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
				tmptri[2].sides[2] = true;
				tmptri[3].sides[1] = true;
				_node[_nodenum].type = true;
			}

			tmptri[0].getangle(_node);
			tmptri[1].getangle(_node);
			tmptri[2].getangle(_node);
			tmptri[3].getangle(_node);

			stack.push_back(_nowtri);
			stack.push_back(neitri);
			stack.push_back(_element.size());
			stack.push_back(_element.size() + 1);

			_element[_nowtri].copy(tmptri[0]);
			_element[neitri].copy(tmptri[1]);
			_element.push_back(tmptri[2]);
			_element.push_back(tmptri[3]);
		}
		//�ӂ�����ŗאڗv�f�����݂��Ȃ��Ƃ�
		else {
			std::vector<Element<T> > tmptri(2);

			tmptri[0].setnode(_nodenum, _element[_nowtri].nodes[nownode], _element[_nowtri].nodes[(nownode + 1) % 3]);
			tmptri[0].setneighbor(_element[_nowtri].neighbors[(nownode + 2) % 3], -1, _element.size());
			tmptri[0].setside(_element[_nowtri].sides[(nownode + 2) % 3], false, false);

			tmptri[1].setnode(_nodenum, _element[_nowtri].nodes[(nownode + 2) % 3], _element[_nowtri].nodes[nownode]);
			tmptri[1].setneighbor(_element[_nowtri].neighbors[(nownode + 1) % 3], _nowtri, -1);
			tmptri[1].setside(_element[_nowtri].sides[(nownode + 1) % 3], false, false);

			int nei1 = _element[_nowtri].neighbors[(nownode + 1) % 3];
			if (nei1 != -1) {
				_element[nei1].neighbors[_element[nei1].oppositenode(_nowtri)] = _element.size();
			}

			if (tmptri[0].nodes[1] == _nodenumm1 || tmptri[0].nodes[1] == _nodenump1) {
				tmptri[1].sides[1] = true;
				tmptri[0].sides[2] = true;
			}
			if (tmptri[0].nodes[2] == _nodenumm1 || tmptri[0].nodes[2] == _nodenump1) {
				tmptri[0].sides[1] = true;
			}
			if (tmptri[1].nodes[1] == _nodenumm1 || tmptri[1].nodes[1] == _nodenump1) {
				tmptri[1].sides[2] = true;
			}

			if (_element[_nowtri].sides[nownode] == true) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
				_node[_nodenum].type = true;
			}

			tmptri[0].getangle(_node);
			tmptri[1].getangle(_node);

			stack.push_back(_nowtri);
			stack.push_back(_element.size());

			_element[_nowtri].copy(tmptri[0]);
			_element.push_back(tmptri[1]);
		}
		//.....�X���b�s���O.....
		swapping(_node, _element, stack, _nodenump1, _nodenumm1);
	}


	template<class T>
	void Delaunay<T>::swapping(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element, std::vector<int> &_stack, int _nodenump1, int _nodenumm1) {
		while (_stack.size() > 0) {
			//�X�^�b�N�����̗v�f�����o��
			int nowstack = _stack[_stack.size() - 1];
			_stack.pop_back();
			//�X�^�b�N������o�����v�f�ɗאڂ���v�f���擾
			int neighbortri = _element[nowstack].neighbors[0];
			//�אڂ���O�p�`�v�f�����݂���Ƃ�
			if (neighbortri >= 0 && _element[neighbortri].active == true) {
				int neighbornode = _element[neighbortri].oppositenode(nowstack);
				T r0 = _node[_element[nowstack].nodes[1]].distance(_node[_element[nowstack].nodes[2]]);
				T r1 = _node[_element[nowstack].nodes[0]].distance(_node[_element[neighbortri].nodes[neighbornode]]);
				if (r0 > r1
					&& _element[nowstack].inouton(_element[neighbortri].nodes[neighbornode], _node) == -1
					&& _element[neighbortri].inouton(_element[nowstack].nodes[0], _node) == -(neighbornode + 1)
					&& _element[nowstack].sides[0] == false) {
					
					Element<T> tmpelement;
					tmpelement.copy(_element[neighbortri]);

					int neighbor1 = tmpelement.neighbors[(neighbornode + 1) % 3];
					if (neighbor1 >= 0) {
						_element[neighbor1].neighbors[_element[neighbor1].oppositenode(neighbortri)] = nowstack;
					}

					int neighbor2 = _element[nowstack].neighbors[1];
					if (neighbor2 >= 0) {
						_element[neighbor2].neighbors[_element[neighbor2].oppositenode(nowstack)] = neighbortri;
					}

					_element[neighbortri].setside(tmpelement.sides[(neighbornode + 2) % 3], _element[nowstack].sides[1], false);
					_element[neighbortri].setnode(_element[nowstack].nodes[0], tmpelement.nodes[neighbornode], _element[nowstack].nodes[2]);
					_element[neighbortri].setneighbor(tmpelement.neighbors[(neighbornode + 2) % 3], _element[nowstack].neighbors[1], nowstack);

					_element[nowstack].setside(tmpelement.sides[(neighbornode + 1) % 3], false, _element[nowstack].sides[2]);
					_element[nowstack].setnode(_element[nowstack].nodes[0], _element[nowstack].nodes[1], tmpelement.nodes[neighbornode]);
					_element[nowstack].setneighbor(tmpelement.neighbors[(neighbornode + 1) % 3], neighbortri, _element[nowstack].neighbors[2]);

					//���E�Ӕ���
					if (_element[nowstack].nodes[2] == _nodenumm1 || _element[nowstack].nodes[2] == _nodenump1) {
						_element[nowstack].sides[1] = true;
						_element[neighbortri].sides[2] = true;
					}

					_element[nowstack].getangle(_node);
					_element[neighbortri].getangle(_node);

					_stack.push_back(nowstack);
					_stack.push_back(neighbortri);
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::getsupertriangle(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements) {
		//Get distance maximam
		T rmax = T();
		Node<T> o;
		for (auto& node : _nodes) {
			T tmpr = node.distance(o);
			if (rmax < tmpr) {
				rmax = tmpr;
			}
		}
		
		//Generate SuperTriangle
		for (int i = 0; i < 3; i++) {
			_nodes.push_back(Node<T>(-2.0*rmax * sin(2.0*M_PI*i / 3.0), 2.0*rmax * cos(2.0*M_PI*i / 3.0)));
		}
		_elements.push_back(Element<T>(_nodes.size() - 3, _nodes.size() - 2, _nodes.size() - 1));
	}


	template<class T>
	void Delaunay<T>::deletesupertriangle(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element) {
		for (int i = _element.size() - 1; i >= 0; i--) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					if (_element[i].nodes[j] == _node.size() - 1 - k) {
						_element[i].active = false;
						break;
					}
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::getboundary(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element, Boundary _boundary) {
		for (int i = 0; i < _boundary.nodelist.size(); i++) {
			//.....�܂��ݒu����Ă��Ȃ��Ƃ�.....
			if (_node[_boundary.nodelist[i]].set == false) {
				_node[_boundary.nodelist[i]].set = true;
				int nowtri = 0;
				if (_element.size() > 0) {
					nowtri = _element.size() - 1;
				}
				_node[_boundary.nodelist[i]].type = true;

				//.....��܎O�p�`�̒T��.....
				for (int j = 0; j < _element.size(); j++) {
					int pos = _element[nowtri].inouton(_boundary.nodelist[i], _node);
					//�v�f�O�ɂ��鎞
					if (pos < 0 || _element[nowtri].active == false) {
						if (_element[nowtri].neighbors[abs(pos) - 1] >= 0) {
							nowtri = _element[nowtri].neighbors[abs(pos) - 1];
						}
						else {
							std::cout << "Out of triangle Error!\n";
						}
					}
					//�v�f���ɂ��鎞
					else if (pos == 0) {
						if (i == 0) {
							getelementin(_node, _element, nowtri, _boundary.nodelist[i + 1], _boundary.nodelist[i], -2);
						}
						else if (i == _boundary.nodelist.size() - 1) {
							getelementin(_node, _element, nowtri, _boundary.nodelist[0], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
						}
						else {
							getelementin(_node, _element, nowtri, _boundary.nodelist[i + 1], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
						}
						break;
					}
					//�ӏ�ɂ��鎞
					else {
						if (i == 0) {
							getelementon(_node, _element, nowtri, pos - 1, _boundary.nodelist[i + 1], _boundary.nodelist[i], -2);
						}
						else if (i == _boundary.nodelist.size() - 1) {
							getelementon(_node, _element, nowtri, pos - 1, _boundary.nodelist[0], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
						}
						else {
							getelementon(_node, _element, nowtri, pos - 1, _boundary.nodelist[i + 1], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
						}
						break;
					}
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::deactivate(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element, Boundary _boundary) {
		for (int i = _element.size() - 1; i >= 0; i--) {
			int nodeorder[3];
			for (int j = 0; j < 3; j++) {
				nodeorder[j] = _boundary.order(_element[i].nodes[j]);
			}

			if (_boundary.type == true) {
				if (_element[i].check == false) {
					if (nodeorder[0] >= 0 && nodeorder[1] >= 0 && nodeorder[2] >= 0) {
						_element[i].check = true;
						if ((nodeorder[0] < nodeorder[1] && nodeorder[1] < nodeorder[2])
							|| (nodeorder[1] < nodeorder[2] && nodeorder[2] < nodeorder[0])
							|| (nodeorder[2] < nodeorder[0] && nodeorder[0] < nodeorder[1])) {
							_element[i].active = false;
						}
					}
				}
			}
			else {
				if (_element[i].check == false) {
					if ((nodeorder[0] < 0 && nodeorder[1] >= 0 && nodeorder[2] >= 0)
						|| (nodeorder[0] >= 0 && nodeorder[1] < 0 && nodeorder[2] >= 0)
						|| (nodeorder[0] >= 0 && nodeorder[1] >= 0 && nodeorder[2] < 0)
						|| (nodeorder[0] < 0 && nodeorder[1] < 0 && nodeorder[2] >= 0)
						|| (nodeorder[0] >= 0 && nodeorder[1] < 0 && nodeorder[2] < 0)
						|| (nodeorder[0] < 0 && nodeorder[1] >= 0 && nodeorder[2] < 0)) {
						_element[i].check = true;
					}
					else if (nodeorder[0] >= 0 && nodeorder[1] >= 0 && nodeorder[2] >= 0) {
						_element[i].check = true;
						_element[i].active = false;
					}
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::sortelement(std::vector<Element<T> > &_element) {
		for (int i = 0; i < _element.size(); i++) {
			for (int j = 0; j < 3; j++) {
				_element[i].neighbors[j] = -1;
				for (int k = 0; k < _element.size(); k++) {
					if (i != k && ((_element[i].nodes[(j + 1) % 3] == _element[k].nodes[0] && _element[i].nodes[(j + 2) % 3] == _element[k].nodes[2])
						|| (_element[i].nodes[(j + 1) % 3] == _element[k].nodes[1] && _element[i].nodes[(j + 2) % 3] == _element[k].nodes[0])
						|| (_element[i].nodes[(j + 1) % 3] == _element[k].nodes[2] && _element[i].nodes[(j + 2) % 3] == _element[k].nodes[1]))) {
						_element[i].neighbors[j] = k;
						break;
					}
				}
			}
		}
	}


	template<class T>
	void Delaunay<T>::getinternalelement(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element, T _maxside) {
		//�ӂ̒����v�f�𕪊�
		for (int i = 0; i < ADDITIONALNODENUM0; i++) {
			//�ǉ�����ߓ_�𐶐�
			T maxside = _node[_element[0].nodes[0]].distance(_node[_element[0].nodes[1]]);
			int maxelement = 0, maxnode = 0;
			for (int j = 0; j < _element.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (maxside < _node[_element[j].nodes[(k + 1) % 3]].distance(_node[_element[j].nodes[(k + 2) % 3]]) && _element[j].active == true) {
						maxside = _node[_element[j].nodes[(k + 1) % 3]].distance(_node[_element[j].nodes[(k + 2) % 3]]);
						maxelement = j;
						maxnode = k;
					}
				}
			}

			if (maxside < _maxside) {
				break;
			}

			_node.push_back(Node<T>(
				(_node[_element[maxelement].nodes[(maxnode + 1) % 3]].x + _node[_element[maxelement].nodes[(maxnode + 2) % 3]].x) / 2.0,
				(_node[_element[maxelement].nodes[(maxnode + 1) % 3]].y + _node[_element[maxelement].nodes[(maxnode + 2) % 3]].y) / 2.0
			));

			getelementon(_node, _element, maxelement, maxnode, -2, _node.size() - 1, -2);
		}
	}


	template<class T>
	void Delaunay<T>::deleteelement(std::vector<Element<T> > &_element) {
		for (int i = _element.size() - 1; i >= 0; i--) {
			if (_element[i].active == false || _element[i].check == false) {
				_element.erase(_element.begin() + i);
			}
		}
	}


	template<class T>
	void Delaunay<T>::laplacian(std::vector<Node<T> > &_node, std::vector<Element<T> > &_element, int _maxnum) {
		std::vector<int> logstack;
		int logstacknum = 100;			//�����ߓ_�΂�����C�����Ȃ��悤�ɒ��߂ɏC���������̂��X�g�b�N
		for (int i = 0; i < _maxnum; i++) {
			//���p�ő�v�f�Ƃ��̗v�f�ԍ���T��
			T angmax = 0.0;
			int elemax = -1, nodemax = -1;
			for (int j = 0; j < _element.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (angmax < fabs(60.0 - _element[j].angles[k]) && _node[_element[j].nodes[k]].type == false) {
						angmax = fabs(60.0 - _element[j].angles[k]);
						elemax = j;
						nodemax = _element[j].nodes[k];
					}
				}
			}
			if (elemax < 0) {
				break;
			}
			
			_node[nodemax].type = true;
			logstack.push_back(nodemax);
			if (logstack.size() > logstacknum) {
				_node[logstack[0]].type = false;
				logstack.erase(logstack.begin());
			}

			//���p�ő�v�f�Ɨאڂ���v�f���X�^�b�N�ɓ����
			std::vector<int> stack;
			int nowtri = elemax;
			do{
				stack.push_back(nowtri);
				nowtri = _element[nowtri].neighbors[(_element[nowtri].nodeorder(nodemax) + 1) % 3];
			}while (nowtri != elemax);

			//�V�������W�����߂�
			T xdash = 0.0, ydash = 0.0;
			for (int j = 0; j < stack.size(); j++) {
				xdash += _node[_element[stack[j]].nodes[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].x 
					+ _node[_element[stack[j]].nodes[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].x;
				ydash += _node[_element[stack[j]].nodes[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].y
					+ _node[_element[stack[j]].nodes[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].y;
			}
			xdash /= 2.0*(T)stack.size();
			ydash /= 2.0*(T)stack.size();

			_node[nodemax].x = xdash;
			_node[nodemax].y = ydash;

			//���p�̍Čv�Z
			for (int j = 0; j < stack.size(); j++) {
				_element[stack[j]].getangle(_node);
			}
		}
	}
}