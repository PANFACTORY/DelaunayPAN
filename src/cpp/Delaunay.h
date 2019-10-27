//*****************************************************************************
//Title		:src/cpp/Delaunay.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <iostream>
#include <vector>
#include <array>
#define _USE_MATH_DEFINES
#include <cmath>


#include "Node.h"
#include "Element.h"
#include "Boundary.h"


namespace DelaunayPAN{
	#define ADDITIONALNODENUM0	100000			
	

	template<class T>
	void delaunaymain(std::vector<Node<T> > &_nodes, std::vector<Element<T> > &_elements, std::vector<Boundary> &_boundaries, T _maxsize, int _laplaciannum) {
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
		//sortelement(_elements);

		/*if (_maxsize > 0) {
			//----------subdivide----------
			getinternalelement(_nodes, _elements, _maxsize);

			//----------Laplacian smoothing----------
			laplacian(_nodes, _elements, _laplaciannum);
		}*/
	}


	template<class T>
	void getelementin(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements, int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1) {
		std::vector<int> stack;
		std::array<Element<T>, 3> tmpelement;

		tmpelement[0].sides[0] = _elements[_nowtri].sides[0];
		tmpelement[1].sides[0] = _elements[_nowtri].sides[1];
		tmpelement[2].sides[0] = _elements[_nowtri].sides[2];

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

		swapping(_nodes, _elements, stack, _nodenump1, _nodenumm1);
	}


	template<class T>
	void getelementon(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements, int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1) {
		std::vector<int> stack;
		int nownode = _pos;
		int neitri = _elements[_nowtri].neighbors[nownode];
		
		if (neitri != -1 && _elements[neitri].active == true) {
			int neinode = _elements[neitri].oppositenode(_nowtri);
			std::array<Element<T>, 4> tmptri;			//0:nowtri	1:neitri

			tmptri[0].setnode(_nodenum, _elements[_nowtri].nodes[nownode], _elements[_nowtri].nodes[(nownode + 1) % 3]);
			tmptri[0].setneighbor(_elements[_nowtri].neighbors[(nownode + 2) % 3], neitri, _elements.size());
			tmptri[0].setside(_elements[_nowtri].sides[(nownode + 2) % 3], false, false);

			tmptri[1].setnode(_nodenum, _elements[_nowtri].nodes[(nownode + 1) % 3], _elements[neitri].nodes[neinode]);
			tmptri[1].setneighbor(_elements[neitri].neighbors[(neinode + 1) % 3], _elements.size() + 1, _nowtri);
			tmptri[1].setside(_elements[neitri].sides[(neinode + 1) % 3], false, false);

			tmptri[2].setnode(_nodenum, _elements[_nowtri].nodes[(nownode + 2) % 3], _elements[_nowtri].nodes[nownode]);
			tmptri[2].setneighbor(_elements[_nowtri].neighbors[(nownode + 1) % 3], _nowtri, _elements.size() + 1);
			tmptri[2].sides[0] = _elements[_nowtri].sides[(nownode + 1) % 3];

			tmptri[3].setnode(_nodenum, _elements[neitri].nodes[neinode], _elements[_nowtri].nodes[(nownode + 2) % 3]);
			tmptri[3].setneighbor(_elements[neitri].neighbors[(neinode + 2) % 3], _elements.size(), neitri);
			tmptri[3].sides[0] = _elements[neitri].sides[(neinode + 2) % 3];

			int nei1 = _elements[_nowtri].neighbors[(nownode + 1) % 3];
			if (nei1 != -1) {
				_elements[nei1].neighbors[_elements[nei1].oppositenode(_nowtri)] = _elements.size();
			}

			int nei2 = _elements[neitri].neighbors[(neinode + 2) % 3];
			if (nei2 != -1) {
				_elements[nei2].neighbors[_elements[nei2].oppositenode(neitri)] = _elements.size() + 1;
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

			if (_elements[_nowtri].sides[nownode] == true) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
				tmptri[2].sides[2] = true;
				tmptri[3].sides[1] = true;
				_nodes[_nodenum].isonboundary = true;
			}

			tmptri[0].getangle(_nodes);
			tmptri[1].getangle(_nodes);
			tmptri[2].getangle(_nodes);
			tmptri[3].getangle(_nodes);

			stack.push_back(_nowtri);
			stack.push_back(neitri);
			stack.push_back(_elements.size());
			stack.push_back(_elements.size() + 1);

			_elements[_nowtri].copy(tmptri[0]);
			_elements[neitri].copy(tmptri[1]);
			_elements.push_back(tmptri[2]);
			_elements.push_back(tmptri[3]);
		}
		
		else {
			std::array<Element<T>, 2> tmptri;

			tmptri[0].setnode(_nodenum, _elements[_nowtri].nodes[nownode], _elements[_nowtri].nodes[(nownode + 1) % 3]);
			tmptri[0].setneighbor(_elements[_nowtri].neighbors[(nownode + 2) % 3], -1, _elements.size());
			tmptri[0].setside(_elements[_nowtri].sides[(nownode + 2) % 3], false, false);

			tmptri[1].setnode(_nodenum, _elements[_nowtri].nodes[(nownode + 2) % 3], _elements[_nowtri].nodes[nownode]);
			tmptri[1].setneighbor(_elements[_nowtri].neighbors[(nownode + 1) % 3], _nowtri, -1);
			tmptri[1].setside(_elements[_nowtri].sides[(nownode + 1) % 3], false, false);

			int nei1 = _elements[_nowtri].neighbors[(nownode + 1) % 3];
			if (nei1 != -1) {
				_elements[nei1].neighbors[_elements[nei1].oppositenode(_nowtri)] = _elements.size();
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

			if (_elements[_nowtri].sides[nownode] == true) {
				tmptri[0].sides[1] = true;
				tmptri[1].sides[2] = true;
				_nodes[_nodenum].isonboundary = true;
			}

			tmptri[0].getangle(_nodes);
			tmptri[1].getangle(_nodes);

			stack.push_back(_nowtri);
			stack.push_back(_elements.size());

			_elements[_nowtri].copy(tmptri[0]);
			_elements.push_back(tmptri[1]);
		}
		
		swapping(_nodes, _elements, stack, _nodenump1, _nodenumm1);
	}


	template<class T>
	void swapping(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements, std::vector<int> &_stack, int _nodenump1, int _nodenumm1) {
		while (_stack.size() > 0) {
			int nowstack = _stack[_stack.size() - 1];
			_stack.pop_back();
			
			int neighbortri = _elements[nowstack].neighbors[0];
			
			if (neighbortri >= 0 && _elements[neighbortri].active == true) {
				int neighbornode = _elements[neighbortri].oppositenode(nowstack);
				T r0 = _nodes[_elements[nowstack].nodes[1]].distance(_nodes[_elements[nowstack].nodes[2]]);
				T r1 = _nodes[_elements[nowstack].nodes[0]].distance(_nodes[_elements[neighbortri].nodes[neighbornode]]);
				if (r0 > r1
					&& _elements[nowstack].inouton(_elements[neighbortri].nodes[neighbornode], _nodes) == -1
					&& _elements[neighbortri].inouton(_elements[nowstack].nodes[0], _nodes) == -(neighbornode + 1)
					&& _elements[nowstack].sides[0] == false) {
					
					Element<T> tmpelement = Element<T>(_elements[neighbortri]);

					int neighbor1 = tmpelement.neighbors[(neighbornode + 1) % 3];
					if (neighbor1 >= 0) {
						_elements[neighbor1].neighbors[_elements[neighbor1].oppositenode(neighbortri)] = nowstack;
					}

					int neighbor2 = _elements[nowstack].neighbors[1];
					if (neighbor2 >= 0) {
						_elements[neighbor2].neighbors[_elements[neighbor2].oppositenode(nowstack)] = neighbortri;
					}

					_elements[neighbortri].setside(tmpelement.sides[(neighbornode + 2) % 3], _elements[nowstack].sides[1], false);
					_elements[neighbortri].setnode(_elements[nowstack].nodes[0], tmpelement.nodes[neighbornode], _elements[nowstack].nodes[2]);
					_elements[neighbortri].setneighbor(tmpelement.neighbors[(neighbornode + 2) % 3], _elements[nowstack].neighbors[1], nowstack);

					_elements[nowstack].setside(tmpelement.sides[(neighbornode + 1) % 3], false, _elements[nowstack].sides[2]);
					_elements[nowstack].setnode(_elements[nowstack].nodes[0], _elements[nowstack].nodes[1], tmpelement.nodes[neighbornode]);
					_elements[nowstack].setneighbor(tmpelement.neighbors[(neighbornode + 1) % 3], neighbortri, _elements[nowstack].neighbors[2]);

					if (_elements[nowstack].nodes[2] == _nodenumm1 || _elements[nowstack].nodes[2] == _nodenump1) {
						_elements[nowstack].sides[1] = true;
						_elements[neighbortri].sides[2] = true;
					}

					_elements[nowstack].getangle(_nodes);
					_elements[neighbortri].getangle(_nodes);

					_stack.push_back(nowstack);
					_stack.push_back(neighbortri);
				}
			}
		}
	}


	template<class T>
	void getsupertriangle(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements) {
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
	void deletesupertriangle(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements) {
		for (int i = _elements.size() - 1; i >= 0; i--) {
			for (const auto& node : _elements[i].nodes) {
				if (node == _nodes.size() - 1 || node == _nodes.size() - 2 || node == _nodes.size() - 3) {
					_elements[i].active = false;
					break;
				}
			}
		}
	}


	template<class T>
	void getboundary(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements, Boundary _boundary) {
		for (int i = 0; i < _boundary.nodelists.size(); i++) {
			//.....Add nodes on boundary into district.....
			if (_nodes[_boundary.nodelists[i]].isset == false) {
				_nodes[_boundary.nodelists[i]].isset = true;
				int nowtri = 0;
				if (_elements.size() > 0) {
					nowtri = _elements.size() - 1;
				}
				_nodes[_boundary.nodelists[i]].isonboundary = true;

				//.....Search element in which node is.....
				for (int j = 0; j < _elements.size(); j++) {
					int pos = _elements[nowtri].inouton(_boundary.nodelists[i], _nodes);
					//if not in or out
					if (pos < 0 || _elements[nowtri].active == false) {
						if (_elements[nowtri].neighbors[abs(pos) - 1] >= 0) {
							nowtri = _elements[nowtri].neighbors[abs(pos) - 1];
						} else {
							std::cout << "Out of triangle Error!\n";
						}
					}
					//if in 
					else if (pos == 0) {
						if (i == 0) {
							getelementin(_nodes, _elements, nowtri, _boundary.nodelists[i + 1], _boundary.nodelists[i], -2);
						} else if (i == _boundary.nodelists.size() - 1) {
							getelementin(_nodes, _elements, nowtri, _boundary.nodelists[0], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						} else {
							getelementin(_nodes, _elements, nowtri, _boundary.nodelists[i + 1], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						}
						break;
					}
					//if on
					else {
						if (i == 0) {
							getelementon(_nodes, _elements, nowtri, pos - 1, _boundary.nodelists[i + 1], _boundary.nodelists[i], -2);
						} else if (i == _boundary.nodelists.size() - 1) {
							getelementon(_nodes, _elements, nowtri, pos - 1, _boundary.nodelists[0], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						} else {
							getelementon(_nodes, _elements, nowtri, pos - 1, _boundary.nodelists[i + 1], _boundary.nodelists[i], _boundary.nodelists[i - 1]);
						}
						break;
					}
				}
			}
		}
	}


	template<class T>
	void deactivate(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements, Boundary _boundary) {
		for (auto& element : _elements) {
			if(element.check == false){
				//.....get order of node on boundary.....
				std::array<int, 3> nodeorders{ _boundary.order(element.nodes[0]), _boundary.order(element.nodes[1]), _boundary.order(element.nodes[2]) };
				
				//.....external boundary.....
				if (_boundary.type == true) {
					if (nodeorders[0] >= 0 && nodeorders[1] >= 0 && nodeorders[2] >= 0) {
						element.check = true;
						if ((nodeorders[0] < nodeorders[1] && nodeorders[1] < nodeorders[2])
							|| (nodeorders[1] < nodeorders[2] && nodeorders[2] < nodeorders[0])
							|| (nodeorders[2] < nodeorders[0] && nodeorders[0] < nodeorders[1])) {
							element.active = false;
						}
					}
				}
				//.....internal boundary.....
				else {
					//.....some nodes are on other boundary.....
					if ((nodeorders[0] < 0 && nodeorders[1] >= 0 && nodeorders[2] >= 0)
						|| (nodeorders[0] >= 0 && nodeorders[1] < 0 && nodeorders[2] >= 0)
						|| (nodeorders[0] >= 0 && nodeorders[1] >= 0 && nodeorders[2] < 0)
						|| (nodeorders[0] < 0 && nodeorders[1] < 0 && nodeorders[2] >= 0)
						|| (nodeorders[0] >= 0 && nodeorders[1] < 0 && nodeorders[2] < 0)
						|| (nodeorders[0] < 0 && nodeorders[1] >= 0 && nodeorders[2] < 0)) {
						element.check = true;
					} else if (nodeorders[0] >= 0 && nodeorders[1] >= 0 && nodeorders[2] >= 0) {
						element.check = true;
						if ((nodeorders[0] < nodeorders[1] && nodeorders[1] < nodeorders[2])
							|| (nodeorders[1] < nodeorders[2] && nodeorders[2] < nodeorders[0])
							|| (nodeorders[2] < nodeorders[0] && nodeorders[0] < nodeorders[1])) {
							element.active = false;
						}
					}
				}
			}
		}
	}


	template<class T>
	void sortelement(std::vector<Element<T> >& _elements) {
		for (int i = 0; i < _elements.size(); i++) {
			for (int j = 0; j < 3; j++) {
				_elements[i].neighbors[j] = -1;
				for (int k = 0; k < _elements.size(); k++) {
					if (i != k && ((_elements[i].nodes[(j + 1) % 3] == _elements[k].nodes[0] && _elements[i].nodes[(j + 2) % 3] == _elements[k].nodes[2])
						|| (_elements[i].nodes[(j + 1) % 3] == _elements[k].nodes[1] && _elements[i].nodes[(j + 2) % 3] == _elements[k].nodes[0])
						|| (_elements[i].nodes[(j + 1) % 3] == _elements[k].nodes[2] && _elements[i].nodes[(j + 2) % 3] == _elements[k].nodes[1]))) {
						_elements[i].neighbors[j] = k;
						break;
					}
				}
			}
		}
	}


	template<class T>
	void getinternalelement(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements, T _maxside) {
		//�ӂ̒����v�f�𕪊�
		for (int i = 0; i < ADDITIONALNODENUM0; i++) {
			//�ǉ�����ߓ_�𐶐�
			T maxside = _nodes[_elements[0].nodes[0]].distance(_nodes[_elements[0].nodes[1]]);
			int maxelement = 0, maxnode = 0;
			for (int j = 0; j < _elements.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (maxside < _nodes[_elements[j].nodes[(k + 1) % 3]].distance(_nodes[_elements[j].nodes[(k + 2) % 3]]) && _elements[j].active == true) {
						maxside = _nodes[_elements[j].nodes[(k + 1) % 3]].distance(_nodes[_elements[j].nodes[(k + 2) % 3]]);
						maxelement = j;
						maxnode = k;
					}
				}
			}

			if (maxside < _maxside) {
				break;
			}

			_nodes.push_back(Node<T>(
				(_nodes[_elements[maxelement].nodes[(maxnode + 1) % 3]].x + _nodes[_elements[maxelement].nodes[(maxnode + 2) % 3]].x) / 2.0,
				(_nodes[_elements[maxelement].nodes[(maxnode + 1) % 3]].y + _nodes[_elements[maxelement].nodes[(maxnode + 2) % 3]].y) / 2.0
			));

			getelementon(_nodes, _elements, maxelement, maxnode, -2, _nodes.size() - 1, -2);
		}
	}


	template<class T>
	void deleteelement(std::vector<Element<T> >& _elements) {
		for (int i = _elements.size() - 1; i >= 0; i--) {
			if (_elements[i].active == false || _elements[i].check == false) {
				_elements.erase(_elements.begin() + i);
			}
		}
	}


	template<class T>
	void laplacian(std::vector<Node<T> >& _nodes, std::vector<Element<T> >& _elements, int _maxnum) {
		std::vector<int> logstack;
		int logstacknum = 100;			
		for (int i = 0; i < _maxnum; i++) {
			T angmax = 0.0;
			int elemax = -1, nodemax = -1;
			for (int j = 0; j < _elements.size(); j++) {
				for (int k = 0; k < 3; k++) {
					if (angmax < fabs(60.0 - _elements[j].angles[k]) && _nodes[_elements[j].nodes[k]].isonboundary == false) {
						angmax = fabs(60.0 - _elements[j].angles[k]);
						elemax = j;
						nodemax = _elements[j].nodes[k];
					}
				}
			}
			if (elemax < 0) {
				break;
			}
			
			_nodes[nodemax].isonboundary = true;
			logstack.push_back(nodemax);
			if (logstack.size() > logstacknum) {
				_nodes[logstack[0]].isonboundary = false;
				logstack.erase(logstack.begin());
			}

			std::vector<int> stack;
			int nowtri = elemax;
			do{
				stack.push_back(nowtri);
				nowtri = _elements[nowtri].neighbors[(_elements[nowtri].nodeorder(nodemax) + 1) % 3];
			}while (nowtri != elemax);

			T xdash = T(), ydash = T();
			for (auto j : stack) {
				xdash += _nodes[_elements[j].nodes[(_elements[j].nodeorder(nodemax) + 1) % 3]].x 
					+ _nodes[_elements[j].nodes[(_elements[j].nodeorder(nodemax) + 1) % 3]].x;
				ydash += _nodes[_elements[j].nodes[(_elements[j].nodeorder(nodemax) + 1) % 3]].y
					+ _nodes[_elements[j].nodes[(_elements[j].nodeorder(nodemax) + 1) % 3]].y;
			}
			xdash /= 2.0*(T)stack.size();
			ydash /= 2.0*(T)stack.size();

			_nodes[nodemax].x = xdash;
			_nodes[nodemax].y = ydash;

			for (auto j : stack) {
				_elements[j].getangle(_nodes);
			}
		}
	}
}