//*****************************************************************************
//Title		:Delaunay SolverClass
//Purpose	:Solver for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/27
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************

#include "pch.h"
#include "DelaunayClass.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>


DelaunayClass::DelaunayClass(){}


DelaunayClass::~DelaunayClass(){}


//*****************************************************************************
//SuperTriangleの生成
//*****************************************************************************
void DelaunayClass::getsupertriangle(vector<NodeClass> &_node, vector<ElementClass> &_element) {
	//原点から最も離れた点を探索
	double rmax = 0.0;
	for (int i = 0; i < _node.size(); i++) {
		if (rmax < sqrt(pow(_node[i].x, 2.0) + pow(_node[i].y, 2.0))) {
			rmax = sqrt(pow(_node[i].x, 2.0) + pow(_node[i].y, 2.0));
		}
	}
	//rmaxの1.5倍の直径を持つ円を内接円に持つ三角形を生成
	vector<NodeClass> st(3);						
	ElementClass superelement;
	for (int i = 0; i < 3; i++) {
		st[i].x = -2.0*rmax * sin(2.0*M_PI*i / 3.0);
		st[i].y = 2.0*rmax * cos(2.0*M_PI*i / 3.0);
		_node.push_back(st[i]);
		superelement.node[i] = _node.size() - 1;
	}
	_element.push_back(superelement);
}


//*****************************************************************************
//SuperTriangleの除去
//*****************************************************************************
void DelaunayClass::deletesupertriangle(vector<NodeClass> &_node, vector<ElementClass> &_element) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				if (_element[i].node[j] == _node.size() - 1 - k) {
					_element.erase(_element.begin() + i);
					break;
				}
			}
		}
	}
}


//*****************************************************************************
//DelaunayTriangulationによる要素生成
//*****************************************************************************
void DelaunayClass::getelement(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary) {
	int nowtri = 0;															//位置関係を調べている要素
	for (int i = 0; i < _boundary.nodelist.size(); i++) {
		//.....既に設置された点か確認.....
		if (_node[_boundary.nodelist[i]].set == true) {
			break;
		}
		_node[_boundary.nodelist[i]].set = true;

		vector<int> stack;													//スワッピング対象の要素配列
		//.....包含三角形の探索.....
		for (int j = 0; j < _element.size(); j++) {
			int pos = _element[nowtri].inouton(_boundary.nodelist[i], _node);						
			//要素外にある時
			if (pos < 0) {
				if (_element[nowtri].neighbor[abs(pos) - 1] >= 0) {
					nowtri = _element[nowtri].neighbor[abs(pos) - 1];		
				}
				else {
					cout << "Out of triangle Error!\n";
				}
			}
			//要素内にある時
			else if (pos == 0) {
				cout << "in->" << nowtri << "\n";
				vector<ElementClass> tmpelement(2);

				tmpelement[0].side[0] = _element[nowtri].side[1];
				tmpelement[1].side[0] = _element[nowtri].side[2];
				_element[nowtri].setside(_element[nowtri].side[0], false, false);

				//境界辺判定
				if (i > 0) {
					if (_element[nowtri].node[0] == _boundary.nodelist[i - 1]
						|| (i == _boundary.nodelist.size() - 1 && _element[nowtri].node[0] == _boundary.nodelist[0])) {
						tmpelement[0].side[1] = true;
						tmpelement[1].side[2] = true;
					}
					if (_element[nowtri].node[1] == _boundary.nodelist[i - 1]
						|| (i == _boundary.nodelist.size() - 1 && _element[nowtri].node[1] == _boundary.nodelist[0])) {
						tmpelement[1].side[1] = true;
						_element[nowtri].side[2] = true;
					}
					if (_element[nowtri].node[2] == _boundary.nodelist[i - 1]
						|| (i == _boundary.nodelist.size() - 1 && _element[nowtri].node[2] == _boundary.nodelist[0])) {
						_element[nowtri].side[1] = true;
						tmpelement[0].side[2] = true;
					}
				}

				tmpelement[0].setnode(_boundary.nodelist[i], _element[nowtri].node[2], _element[nowtri].node[0]);
				tmpelement[0].setneighbor(_element[nowtri].neighbor[1], _element.size() + 1, nowtri);

				tmpelement[1].setnode(_boundary.nodelist[i], _element[nowtri].node[0], _element[nowtri].node[1]);
				tmpelement[1].setneighbor(_element[nowtri].neighbor[2], nowtri, _element.size());

				for (int k = 0; k < 2; k++) {
					int neighbor = _element[nowtri].neighbor[1 + k];
					if (neighbor >= 0) {
						_element[neighbor].neighbor[_element[neighbor].oppositenode(nowtri)] = _element.size() + k;
					}
				}

				_element[nowtri].setnode(_boundary.nodelist[i], _element[nowtri].node[1], _element[nowtri].node[2]);
				_element[nowtri].setneighbor(_element[nowtri].neighbor[0], _element.size(), _element.size() + 1);

				stack.push_back(nowtri);
				stack.push_back(_element.size());
				stack.push_back(_element.size() + 1);

				_element.push_back(tmpelement[0]);
				_element.push_back(tmpelement[1]);
				break;
			}
			//辺上にある時
			//*******************要修正********************
			else {
				/*int neighborelement = element[nowtri].neighbor[pos - 1];	//辺を挟む要素

				if (neighborelement >= 0) {
					int neighbornode = element[neighborelement].oppositenode(nowtri);	//辺を挟む頂点

					vector<ElementClass> tmpelement(2);
					tmpelement[0].setnode(i, element[neighborelement].node[neighbornode], element[nowtri].node[(pos + 1) % 3]);
					tmpelement[0].setneighbor(element[neighborelement].neighbor[(neighbornode + 2) % 3], element.size() + 1, neighborelement);

					tmpelement[1].setnode(i, element[nowtri].node[(pos + 1) % 3], element[nowtri].node[pos - 1]);
					tmpelement[1].setneighbor(element[nowtri].neighbor[pos % 3], nowtri, element.size());

					for (int k = 0; k < 3; k++) {
						if (element[element[nowtri].neighbor[(pos + 1) % 3]].neighbor[k] == nowtri) {
							element[element[nowtri].neighbor[(pos + 1) % 3]].neighbor[k] = element.size() + 1;
						}
						if (element[element[neighborelement].neighbor[(neighbornode + 2) % 3]].neighbor[k] == neighborelement) {
							element[element[neighborelement].neighbor[(neighbornode + 2) % 3]].neighbor[k] = element.size();
						}
					}

					element[nowtri].setnode(i, element[nowtri].node[pos - 1], element[nowtri].node[pos % 3]);
					element[nowtri].setneighbor(element[nowtri].neighbor[(pos + 1) % 3], neighborelement, element.size() + 1);

					element[neighborelement].setnode(i, element[nowtri].node[pos % 3], element[neighborelement].node[neighbornode]);
					element[neighborelement].setneighbor(element[neighborelement].neighbor[(neighbornode + 1) % 3], neighborelement, element.size());

					element.push_back(tmpelement[0]);
					element.push_back(tmpelement[1]);
				}
				break;*/
			}
		}

		//.....スワッピング.....
		while (stack.size() > 0) {
			//スタック末尾の要素を取り出す
			int nowstack = stack[stack.size() - 1];
			stack.pop_back();
			//スタックから取り出した要素に隣接する要素を取得
			int neighbortri = _element[nowstack].neighbor[0];
			//隣接する三角形要素が存在するとき
			if (neighbortri >= 0) {
				int neighbornode = _element[neighbortri].oppositenode(nowstack);
				double r0 = _node[_element[nowstack].node[1]].distance(_node[_element[nowstack].node[2]]);
				double r1 = _node[_element[nowstack].node[0]].distance(_node[_element[neighbortri].node[neighbornode]]);
				if (r0 > r1
					&& _element[nowstack].inouton(_element[neighbortri].node[neighbornode], _node) == -1
					&& _element[neighbortri].inouton(_element[nowstack].node[0], _node) == -(neighbornode + 1)
					&& _element[nowstack].side[0] == false) {

					ElementClass tmpelement;
					tmpelement.copy(_element[neighbortri]);

					int neighbor1 = tmpelement.neighbor[(neighbornode + 1) % 3];
					if (neighbor1 >= 0) {
						_element[neighbor1].neighbor[_element[neighbor1].oppositenode(neighbortri)] = nowstack;
					}

					int neighbor2 = _element[nowstack].neighbor[1];
					if (neighbor2 >= 0) {
						_element[neighbor2].neighbor[_element[neighbor2].oppositenode(nowstack)] = neighbortri;
					}

					_element[neighbortri].setside(tmpelement.side[(neighbornode + 2) % 3], _element[nowstack].side[1], false);
					_element[neighbortri].setnode(_element[nowstack].node[0], tmpelement.node[neighbornode], _element[nowstack].node[2]);
					_element[neighbortri].setneighbor(tmpelement.neighbor[(neighbornode + 2) % 3], _element[nowstack].neighbor[1], nowstack);

					_element[nowstack].setside(tmpelement.side[(neighbornode + 1) % 3], false, _element[nowstack].side[2]);
					_element[nowstack].setnode(_element[nowstack].node[0], _element[nowstack].node[1], tmpelement.node[neighbornode]);
					_element[nowstack].setneighbor(tmpelement.neighbor[(neighbornode + 1) % 3], neighbortri, _element[nowstack].neighbor[2]);

					//境界辺判定
					if (i > 0) {
						if (_element[nowstack].node[2] == _boundary.nodelist[i - 1]
							|| (i == _boundary.nodelist.size() - 1 && _element[nowstack].node[2] == _boundary.nodelist[0])) {
							_element[nowstack].side[1] = true;
							_element[neighbortri].side[2] = true;
						}
					}

					stack.push_back(nowstack);
					stack.push_back(neighbortri);
				}
			}
		}
	}
}


//*****************************************************************************
//境界外部の要素を削除
//*****************************************************************************
void DelaunayClass::deleteelement(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		int nodeorder[3];
		for (int j = 0; j < 3; j++) {
			nodeorder[j] = _boundary.order(_element[i].node[j]);
		}
		
		if (((nodeorder[0] < nodeorder[1] && nodeorder[1] < nodeorder[2])
			|| (nodeorder[1] < nodeorder[2] && nodeorder[2] < nodeorder[0])
			|| (nodeorder[2] < nodeorder[0] && nodeorder[0] < nodeorder[1]))
			&& (nodeorder[0] >= 0 && nodeorder[1] >= 0 && nodeorder[2] >= 0)) {
			_element.erase(_element.begin() + i);
		}
	}
}