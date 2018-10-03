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


#define ADDITIONALNODENUM0	1000			//細分割数
#define ADDITIONALNODENUM1	0				//細分割数
#define ADDITIONALNODENUM2	0				//細分割数


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
//要素内に点がきたとき
//*****************************************************************************
void DelaunayClass::getelementin(vector<NodeClass> &_node, vector<ElementClass> &_element, int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1) {
	vector<int> stack;
	
	vector<ElementClass> tmpelement(2);

	tmpelement[0].side[0] = _element[_nowtri].side[1];
	tmpelement[1].side[0] = _element[_nowtri].side[2];
	_element[_nowtri].setside(_element[_nowtri].side[0], false, false);

	//境界辺判定
	if (_element[_nowtri].node[0] == _nodenumm1 || _element[_nowtri].node[0] == _nodenump1) {
		tmpelement[0].side[1] = true;
		tmpelement[1].side[2] = true;
	}
	if (_element[_nowtri].node[1] == _nodenumm1 || _element[_nowtri].node[1] == _nodenump1) {
		tmpelement[1].side[1] = true;
		_element[_nowtri].side[2] = true;
	}
	if (_element[_nowtri].node[2] == _nodenumm1 || _element[_nowtri].node[2] == _nodenump1) {
		_element[_nowtri].side[1] = true;
		tmpelement[0].side[2] = true;
	}

	tmpelement[0].setnode(_nodenum, _element[_nowtri].node[2], _element[_nowtri].node[0]);
	tmpelement[0].setneighbor(_element[_nowtri].neighbor[1], _element.size() + 1, _nowtri);

	tmpelement[1].setnode(_nodenum, _element[_nowtri].node[0], _element[_nowtri].node[1]);
	tmpelement[1].setneighbor(_element[_nowtri].neighbor[2], _nowtri, _element.size());

	for (int k = 0; k < 2; k++) {
		int neighbor = _element[_nowtri].neighbor[1 + k];
		if (neighbor >= 0) {
			_element[neighbor].neighbor[_element[neighbor].oppositenode(_nowtri)] = _element.size() + k;
		}
	}

	_element[_nowtri].setnode(_nodenum, _element[_nowtri].node[1], _element[_nowtri].node[2]);
	_element[_nowtri].setneighbor(_element[_nowtri].neighbor[0], _element.size(), _element.size() + 1);

	_element[_nowtri].getangle(_node);
	tmpelement[0].getangle(_node);
	tmpelement[1].getangle(_node);

	stack.push_back(_nowtri);
	stack.push_back(_element.size());
	stack.push_back(_element.size() + 1);

	_element.push_back(tmpelement[0]);
	_element.push_back(tmpelement[1]);

	//.....スワッピング.....
	swapping(_node, _element, stack, _nodenump1, _nodenumm1);
}


//*****************************************************************************
//辺上に点がきたとき
//*****************************************************************************
void DelaunayClass::getelementon(vector<NodeClass> &_node, vector<ElementClass> &_element, int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1) {
	vector<int> stack;
	int nownode = _pos - 1;
	int neitri = _element[_nowtri].neighbor[nownode];
	//辺を挟んで隣接要素が存在するとき
	if (neitri != -1 && _element[neitri].active == true) {
		int neinode = _element[neitri].oppositenode(_nowtri);
		vector<ElementClass> tmptri(4);
		tmptri[0].copy(_element[_nowtri]);
		tmptri[1].copy(_element[neitri]);

		_element[_nowtri].setnode(_nodenum, tmptri[0].node[nownode], tmptri[0].node[(nownode + 1) % 3]);
		_element[_nowtri].setneighbor(tmptri[0].neighbor[(nownode + 2) % 3], neitri, _element.size());
		_element[_nowtri].setside(tmptri[0].side[(nownode + 2) % 3], false, false);

		_element[neitri].setnode(_nodenum, tmptri[0].node[(nownode + 1) % 3], tmptri[1].node[neinode]);
		_element[neitri].setneighbor(tmptri[1].neighbor[(neinode + 1) % 3], _element.size() + 1, _nowtri);
		_element[neitri].setside(tmptri[1].side[(neinode + 1) % 3], false, false);

		tmptri[2].setnode(_nodenum, tmptri[0].node[(nownode + 2) % 3], tmptri[0].node[nownode]);
		tmptri[2].setneighbor(tmptri[0].neighbor[(nownode + 1) % 3], _nowtri, _element.size() + 1);
		tmptri[2].side[0] = tmptri[0].side[(nownode + 1) % 3];

		tmptri[3].setnode(_nodenum, tmptri[1].node[neinode], tmptri[0].node[(nownode + 2) % 3]);
		tmptri[3].setneighbor(tmptri[1].neighbor[(neinode + 2) % 3], _element.size(), neitri);
		tmptri[3].side[0] = tmptri[1].side[(neinode + 2) % 3];

		int nei1 = tmptri[0].neighbor[(nownode + 1) % 3];
		if (nei1 != -1) {
			_element[nei1].neighbor[_element[nei1].oppositenode(_nowtri)] = _element.size();
		}

		int nei2 = tmptri[1].neighbor[(neinode + 2) % 3];
		if (nei2 != -1) {
			_element[nei2].neighbor[_element[nei2].oppositenode(neitri)] = _element.size() + 1;
		}

		if (_element[_nowtri].node[1] == _nodenumm1 || _element[_nowtri].node[1] == _nodenump1) {
			tmptri[2].side[1] = true;
			_element[_nowtri].side[2] = true;
		}
		if (_element[neitri].node[1] == _nodenumm1 || _element[neitri].node[1] == _nodenump1) {
			_element[_nowtri].side[1] = true;
			_element[neitri].side[2] = true;
		}
		if (tmptri[2].node[1] == _nodenumm1 || tmptri[2].node[1] == _nodenump1) {
			tmptri[2].side[2] = true;
			tmptri[3].side[1] = true;
		}
		if (tmptri[3].node[1] == _nodenumm1 || tmptri[2].node[1] == _nodenump1) {
			tmptri[3].side[2] = true;
			_element[neitri].side[1] = true;
		}

		if (tmptri[0].side[nownode] == true) {
			_element[_nowtri].side[1] = true;
			_element[neitri].side[2] = true;
			tmptri[2].side[2] = true;
			tmptri[3].side[1] = true;
		}

		_element[_nowtri].getangle(_node);
		_element[neitri].getangle(_node);
		tmptri[2].getangle(_node);
		tmptri[3].getangle(_node);

		stack.push_back(_nowtri);
		stack.push_back(neitri);
		stack.push_back(_element.size());
		stack.push_back(_element.size() + 1);

		_element.push_back(tmptri[2]);
		_element.push_back(tmptri[3]);
	}
	//辺を挟んで隣接要素が存在しないとき
	else {
		vector<ElementClass> tmptri(2);
		tmptri[0].copy(_element[_nowtri]);

		_element[_nowtri].setnode(_nodenum, tmptri[0].node[nownode], tmptri[0].node[(nownode + 1) % 3]);
		_element[_nowtri].setneighbor(tmptri[0].neighbor[(nownode + 2) % 3], -1, _element.size());
		_element[_nowtri].setside(tmptri[0].side[(nownode + 2) % 3], false, false);

		tmptri[1].setnode(_nodenum, tmptri[0].node[(nownode + 2) % 3], tmptri[0].node[nownode]);
		tmptri[1].setneighbor(tmptri[0].neighbor[(nownode + 1) % 3], _nowtri, -1);
		tmptri[1].setside(tmptri[0].side[(nownode + 1) % 3], false, false);

		int nei1 = tmptri[0].neighbor[(nownode + 1) % 3];
		if (nei1 != -1) {
			_element[nei1].neighbor[_element[nei1].oppositenode(_nowtri)] = _element.size();
		}

		if (_element[_nowtri].node[1] == _nodenumm1 || _element[_nowtri].node[1] == _nodenump1) {
			tmptri[1].side[1] = true;
			_element[_nowtri].side[2] = true;
		}
		if (_element[_nowtri].node[2] == _nodenumm1 || _element[_nowtri].node[2] == _nodenump1) {
			_element[_nowtri].side[1] = true;
		}
		if (tmptri[1].node[1] == _nodenumm1 || tmptri[1].node[1] == _nodenump1) {
			tmptri[1].side[2] = true;
		}

		if (tmptri[0].side[nownode] == true) {
			_element[_nowtri].side[1] = true;
			tmptri[1].side[2] = true;
		}

		_element[_nowtri].getangle(_node);
		tmptri[1].getangle(_node);

		stack.push_back(_nowtri);
		stack.push_back(_element.size());

		_element.push_back(tmptri[1]);
	}
	//.....スワッピング.....
	swapping(_node, _element, stack, _nodenump1, _nodenumm1);
}


//*****************************************************************************
//スワッピング
//*****************************************************************************
void DelaunayClass::swapping(vector<NodeClass> &_node, vector<ElementClass> &_element, vector<int> &_stack, int _nodenump1, int _nodenumm1) {
	while (_stack.size() > 0) {
		//スタック末尾の要素を取り出す
		int nowstack = _stack[_stack.size() - 1];
		_stack.pop_back();
		//スタックから取り出した要素に隣接する要素を取得
		int neighbortri = _element[nowstack].neighbor[0];
		//隣接する三角形要素が存在するとき
		if (neighbortri >= 0 && _element[neighbortri].active == true) {
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
				if (_element[nowstack].node[2] == _nodenumm1 || _element[nowstack].node[2] == _nodenump1) {
					_element[nowstack].side[1] = true;
					_element[neighbortri].side[2] = true;
				}

				_element[nowstack].getangle(_node);
				_element[neighbortri].getangle(_node);

				_stack.push_back(nowstack);
				_stack.push_back(neighbortri);
			}
		}
	}
}


//*****************************************************************************
//境界の生成
//*****************************************************************************
void DelaunayClass::getboundary(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary) {
	for (int i = 0; i < _boundary.nodelist.size(); i++) {
		int nowtri = 0;
		if (_element.size() > 0) {
			nowtri = _element.size() - 1;
		}

		//.....包含三角形の探索.....
		for (int j = 0; j < _element.size(); j++) {
			int pos = _element[nowtri].inouton(_boundary.nodelist[i], _node);
			//要素外にある時
			if (pos < 0 || _element[nowtri].active == false) {
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
				if (i == 0) {
					getelementin(_node, _element, nowtri, _boundary.nodelist[i + 1], _boundary.nodelist[i], -2);
				}
				else if(i == _boundary.nodelist.size() - 1){
					getelementin(_node, _element, nowtri, _boundary.nodelist[0], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
				}
				else {
					getelementin(_node, _element, nowtri, _boundary.nodelist[i + 1], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
				}
				break;
			}
			//辺上にある時
			else {
				cout << "on->" << nowtri << "\n";
				if (i == 0) {
					getelementon(_node, _element, nowtri, pos, _boundary.nodelist[i + 1], _boundary.nodelist[i], -2);
				}
				else if (i == _boundary.nodelist.size() - 1) {
					getelementon(_node, _element, nowtri, pos, _boundary.nodelist[0], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
				}
				else {
					getelementon(_node, _element, nowtri, pos, _boundary.nodelist[i + 1], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
				}
				break;
			}
		}
	}
}


//*****************************************************************************
//境界外部の要素を無効化
//*****************************************************************************
void DelaunayClass::deactivate(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		int nodeorder[3];
		for (int j = 0; j < 3; j++) {
			nodeorder[j] = _boundary.order(_element[i].node[j]);
		}
		
		if (((nodeorder[0] < nodeorder[1] && nodeorder[1] < nodeorder[2])
			|| (nodeorder[1] < nodeorder[2] && nodeorder[2] < nodeorder[0])
			|| (nodeorder[2] < nodeorder[0] && nodeorder[0] < nodeorder[1]))
			&& (nodeorder[0] >= 0 && nodeorder[1] >= 0 && nodeorder[2] >= 0)) {
			_element[i].active = false;
		}
	}
}


//*****************************************************************************
//要素要素間隣接関係の修復
//*****************************************************************************
void DelaunayClass::sortelement(vector<ElementClass> &_element) {
	for (int i = 0; i < _element.size(); i++) {
		for (int j = 0; j < 3; j++) {
			_element[i].neighbor[j] = -1;
			for (int k = 0; k < _element.size(); k++) {
				if (i != k && ((_element[i].node[(j + 1) % 3] == _element[k].node[0] && _element[i].node[(j + 2) % 3] == _element[k].node[2])
					|| (_element[i].node[(j + 1) % 3] == _element[k].node[1] && _element[i].node[(j + 2) % 3] == _element[k].node[0])
					|| (_element[i].node[(j + 1) % 3] == _element[k].node[2] && _element[i].node[(j + 2) % 3] == _element[k].node[1]))) {
					_element[i].neighbor[j] = k;
					break;
				}
			}
		}
	}
}


//*****************************************************************************
//内部点の生成
//*****************************************************************************
void DelaunayClass::getinternalelement(vector<NodeClass> &_node, vector<ElementClass> &_element) {
	//辺の長い要素を分割
	for (int i = 0; i < ADDITIONALNODENUM0; i++) {
		//追加する節点を生成
		double maxside = _node[_element[0].node[0]].distance(_node[_element[0].node[1]]);
		int maxelement = 0, maxnode = 0;
		for (int j = 0; j < _element.size(); j++) {
			for (int k = 0; k < 3; k++) {
				if (maxside < _node[_element[j].node[(k + 1) % 3]].distance(_node[_element[j].node[(k + 2) % 3]]) && _element[j].active == true) {
					maxside = _node[_element[j].node[(k + 1) % 3]].distance(_node[_element[j].node[(k + 2) % 3]]);
					maxelement = j;
					maxnode = k;
				}
			}
		}

		NodeClass addnode;
		addnode.x = (_node[_element[maxelement].node[(maxnode + 1) % 3]].x + _node[_element[maxelement].node[(maxnode + 2) % 3]].x) / 2.0;
		addnode.y = (_node[_element[maxelement].node[(maxnode + 1) % 3]].y + _node[_element[maxelement].node[(maxnode + 2) % 3]].y) / 2.0;
		_node.push_back(addnode);

		//追加した節点で要素分割
		cout << "on-->" << maxelement << "\n";
		getelementon(_node, _element, maxelement, maxnode + 1, -2, _node.size() - 1, -2);
	}
	
	//面積の大きい要素を重心で分割
	for (int i = 0; i < ADDITIONALNODENUM1; i++) {
		double maxspace = _element[0].space(_node);
		int maxelement = 0;
		for (int j = 1; j < _element.size(); j++) {
			if (maxspace < _element[j].space(_node) && _element[j].active == true) {
				maxspace = _element[j].space(_node);
				maxelement = j;
			}
		}
		NodeClass addnode;
		addnode.x = (_node[_element[maxelement].node[0]].x + _node[_element[maxelement].node[1]].x + _node[_element[maxelement].node[2]].x) / 3.0;
		addnode.y = (_node[_element[maxelement].node[0]].y + _node[_element[maxelement].node[1]].y + _node[_element[maxelement].node[2]].y) / 3.0;
		_node.push_back(addnode);

		//追加した節点で要素分割
		cout << "in-->" << maxelement << "\n";
		getelementin(_node, _element, maxelement, -2, _node.size() - 1, -2);
	}
	//ひずみの大きい要素を分割
	for (int i = 0; i < ADDITIONALNODENUM2; i++) {
		//追加する節点を生成
		double maxangle = _element[0].angle[0];
		int maxelement = 0, maxnode = 0;
		for (int j = 0; j < _element.size(); j++) {
			for (int k = 0; k < 3; k++) {
				if (maxangle > _element[j].angle[k] && _element[j].active == true) {
					maxangle = _element[j].angle[k];
					maxelement = j;
					maxnode = k;
				}
			}
		}
		
		NodeClass addnode;
		double a = _node[_element[maxelement].node[maxnode]].distance(_node[_element[maxelement].node[(maxnode + 1) % 3]]);
		double b = _node[_element[maxelement].node[maxnode]].distance(_node[_element[maxelement].node[(maxnode + 2) % 3]]);
		addnode.x = (b * _node[_element[maxelement].node[(maxnode + 1) % 3]].x + a * _node[_element[maxelement].node[(maxnode + 2) % 3]].x) / (a + b);
		addnode.y = (b * _node[_element[maxelement].node[(maxnode + 1) % 3]].y + a * _node[_element[maxelement].node[(maxnode + 2) % 3]].y) / (a + b);
		_node.push_back(addnode);

		//追加した節点で要素分割
		cout << "on-->" << maxelement << "\n";
		getelementon(_node, _element, maxelement, maxnode + 1, -2, _node.size() - 1, -2);
	}
}


//*****************************************************************************
//無効化された要素を削除
//*****************************************************************************
void DelaunayClass::deleteelement(vector<ElementClass> &_element) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		if (_element[i].active == false) {
			_element.erase(_element.begin() + i);
		}
	}
}