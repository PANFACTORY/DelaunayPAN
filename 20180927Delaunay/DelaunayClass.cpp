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

#define ADDITIONALNODENUM1	100				//ç◊ï™äÑêî
#define ADDITIONALNODENUM2	8					//ç◊ï™äÑêî


DelaunayClass::DelaunayClass(){}


DelaunayClass::~DelaunayClass(){}


//*****************************************************************************
//SuperTriangleÇÃê∂ê¨
//*****************************************************************************
void DelaunayClass::getsupertriangle(vector<NodeClass> &_node, vector<ElementClass> &_element) {
	//å¥ì_Ç©ÇÁç≈Ç‡ó£ÇÍÇΩì_ÇíTçı
	double rmax = 0.0;
	for (int i = 0; i < _node.size(); i++) {
		if (rmax < sqrt(pow(_node[i].x, 2.0) + pow(_node[i].y, 2.0))) {
			rmax = sqrt(pow(_node[i].x, 2.0) + pow(_node[i].y, 2.0));
		}
	}
	//rmaxÇÃ1.5î{ÇÃíºåaÇéùÇ¬â~Çì‡ê⁄â~Ç…éùÇ¬éOäpå`Çê∂ê¨
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
//SuperTriangleÇÃèúãé
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
//DelaunayTriangulationÇ…ÇÊÇÈóvëfê∂ê¨
//*****************************************************************************
void DelaunayClass::getelement(vector<NodeClass> &_node, vector<ElementClass> &_element, int _nodenump1, int _nodenum, int _nodenumm1) {
	int nowtri = 0;
	if (_element.size() > 0) {
		nowtri = _element.size() - 1;
	}
	
	vector<int> stack;													//ÉXÉèÉbÉsÉìÉOëŒè€ÇÃóvëfîzóÒ
	//.....ïÔä‹éOäpå`ÇÃíTçı.....
	for (int j = 0; j < _element.size(); j++) {
		int pos = _element[nowtri].inouton(_nodenum, _node);
		//óvëfäOÇ…Ç†ÇÈéû
		if (pos < 0 || _element[nowtri].active == false) {
			if (_element[nowtri].neighbor[abs(pos) - 1] >= 0) {
				nowtri = _element[nowtri].neighbor[abs(pos) - 1];
			}
			else {
				cout << "Out of triangle Error!\n";
			}
		}
		//óvëfì‡Ç…Ç†ÇÈéû
		else if (pos == 0) {
			cout << "in->" << nowtri << "\n";
			vector<ElementClass> tmpelement(2);

			tmpelement[0].side[0] = _element[nowtri].side[1];
			tmpelement[1].side[0] = _element[nowtri].side[2];
			_element[nowtri].setside(_element[nowtri].side[0], false, false);

			//ã´äEï”îªíË
			if (_element[nowtri].node[0] == _nodenumm1 || _element[nowtri].node[0] == _nodenump1) {
				tmpelement[0].side[1] = true;
				tmpelement[1].side[2] = true;
			}
			if (_element[nowtri].node[1] == _nodenumm1 || _element[nowtri].node[1] == _nodenump1) {
				tmpelement[1].side[1] = true;
				_element[nowtri].side[2] = true;
			}
			if (_element[nowtri].node[2] == _nodenumm1 || _element[nowtri].node[2] == _nodenump1) {
				_element[nowtri].side[1] = true;
				tmpelement[0].side[2] = true;
			}

			tmpelement[0].setnode(_nodenum, _element[nowtri].node[2], _element[nowtri].node[0]);
			tmpelement[0].setneighbor(_element[nowtri].neighbor[1], _element.size() + 1, nowtri);

			tmpelement[1].setnode(_nodenum, _element[nowtri].node[0], _element[nowtri].node[1]);
			tmpelement[1].setneighbor(_element[nowtri].neighbor[2], nowtri, _element.size());

			for (int k = 0; k < 2; k++) {
				int neighbor = _element[nowtri].neighbor[1 + k];
				if (neighbor >= 0) {
					_element[neighbor].neighbor[_element[neighbor].oppositenode(nowtri)] = _element.size() + k;
				}
			}

			_element[nowtri].setnode(_nodenum, _element[nowtri].node[1], _element[nowtri].node[2]);
			_element[nowtri].setneighbor(_element[nowtri].neighbor[0], _element.size(), _element.size() + 1);

			_element[nowtri].getangle(_node);
			tmpelement[0].getangle(_node);
			tmpelement[1].getangle(_node);

			stack.push_back(nowtri);
			stack.push_back(_element.size());
			stack.push_back(_element.size() + 1);

			_element.push_back(tmpelement[0]);
			_element.push_back(tmpelement[1]);
			break;
		}
		//ï”è„Ç…Ç†ÇÈéû
		else {
			cout << "on->" << nowtri << "\n";
			int nownode = pos - 1;
			int neitri = _element[nowtri].neighbor[nownode];
			//ï”Çã≤ÇÒÇ≈ó◊ê⁄óvëfÇ™ë∂ç›Ç∑ÇÈÇ∆Ç´
			if (neitri != -1 && _element[neitri].active == true) {
				int neinode = _element[neitri].oppositenode(nowtri);
				vector<ElementClass> tmptri(4);
				tmptri[0].copy(_element[nowtri]);
				tmptri[1].copy(_element[neitri]);

				_element[nowtri].setnode(_nodenum, tmptri[0].node[nownode], tmptri[0].node[(nownode + 1) % 3]);
				_element[nowtri].setneighbor(tmptri[0].neighbor[(nownode + 2) % 3], neitri, _element.size());
				_element[nowtri].setside(tmptri[0].side[(nownode + 2) % 3], false, false);

				_element[neitri].setnode(_nodenum, tmptri[0].node[(nownode + 1) % 3], tmptri[1].node[neinode]);
				_element[neitri].setneighbor(tmptri[1].neighbor[(neinode + 1) % 3], _element.size() + 1, nowtri);
				_element[neitri].setside(tmptri[1].side[(neinode + 1) % 3], false, false);

				tmptri[2].setnode(_nodenum, tmptri[0].node[(nownode + 2) % 3], tmptri[0].node[nownode]);
				tmptri[2].setneighbor(tmptri[0].neighbor[(nownode + 1) % 3], nowtri, _element.size() + 1);
				tmptri[2].side[0] = tmptri[0].side[(nownode + 1) % 3];

				tmptri[3].setnode(_nodenum, tmptri[1].node[neinode], tmptri[0].node[(nownode + 2) % 3]);
				tmptri[3].setneighbor(tmptri[1].neighbor[(neinode + 2) % 3], _element.size(), neitri);
				tmptri[3].side[0] = tmptri[1].side[(neinode + 2) % 3];

				int nei1 = tmptri[0].neighbor[(nownode + 1) % 3];
				if (nei1 != -1) {
					_element[nei1].neighbor[_element[nei1].oppositenode(nowtri)] = _element.size();
				}

				int nei2 = tmptri[1].neighbor[(neinode + 2) % 3];
				if (nei2 != -1) {
					_element[nei2].neighbor[_element[nei2].oppositenode(neitri)] = _element.size() + 1;
				}

				if (_element[nowtri].node[1] == _nodenumm1 || _element[nowtri].node[1] == _nodenump1) {
					tmptri[2].side[1] = true;
					_element[nowtri].side[2] = true;
				}
				if (_element[neitri].node[1] == _nodenumm1 || _element[neitri].node[1] == _nodenump1) {
					_element[nowtri].side[1] = true;
					_element[neitri].side[2] = true;
				}
				if (tmptri[2].node[1] == _nodenumm1 || tmptri[2].node[1] == _nodenump1) {
					tmptri[2].side[2] = true;
					tmptri[3].side[1] = true;
				}
				if (tmptri[3].node[1] == _nodenumm1	|| tmptri[2].node[1] == _nodenump1) {
					tmptri[3].side[2] = true;
					_element[neitri].side[1] = true;
				}

				if (tmptri[0].side[nownode] == true) {
					_element[nowtri].side[1] = true;
					_element[neitri].side[2] = true;
					tmptri[2].side[2] = true;
					tmptri[3].side[1] = true;
				}

				_element[nowtri].getangle(_node);
				_element[neitri].getangle(_node);
				tmptri[2].getangle(_node);
				tmptri[3].getangle(_node);

				_element.push_back(tmptri[2]);
				_element.push_back(tmptri[3]);
			}
			//ï”Çã≤ÇÒÇ≈ó◊ê⁄óvëfÇ™ë∂ç›ÇµÇ»Ç¢Ç∆Ç´
			else {
				vector<ElementClass> tmptri(2);
				tmptri[0].copy(_element[nowtri]);

				_element[nowtri].setnode(_nodenum, tmptri[0].node[nownode], tmptri[0].node[(nownode + 1) % 3]);
				_element[nowtri].setneighbor(tmptri[0].neighbor[(nownode + 2) % 3], -1, _element.size());
				_element[nowtri].setside(tmptri[0].side[(nownode + 2) % 3], false, false);

				tmptri[1].setnode(_nodenum, tmptri[0].node[(nownode + 2) % 3], tmptri[0].node[nownode]);
				tmptri[1].setneighbor(tmptri[0].neighbor[(nownode + 1) % 3], nowtri, -1);
				tmptri[1].setside(tmptri[0].side[(nownode + 1) % 3], false, false);

				int nei1 = tmptri[0].neighbor[(nownode + 1) % 3];
				if (nei1 != -1) {
					_element[nei1].neighbor[_element[nei1].oppositenode(nowtri)] = _element.size();
				}

				if (_element[nowtri].node[1] == _nodenumm1 || _element[nowtri].node[1] == _nodenump1) {
					tmptri[1].side[1] = true;
					_element[nowtri].side[2] = true;
				}
				if (_element[nowtri].node[2] == _nodenumm1 || _element[nowtri].node[2] == _nodenump1) {
					_element[nowtri].side[1] = true;
				}
				if (tmptri[1].node[1] == _nodenumm1	|| tmptri[1].node[1] == _nodenump1) {
					tmptri[1].side[2] = true;
				}

				if (tmptri[0].side[nownode] == true) {
					_element[nowtri].side[1] = true;
					tmptri[1].side[2] = true;
				}

				_element[nowtri].getangle(_node);
				tmptri[1].getangle(_node);

				_element.push_back(tmptri[1]);
			}
			break;
		}
	}

	//.....ÉXÉèÉbÉsÉìÉO.....
	while (stack.size() > 0) {
		//ÉXÉ^ÉbÉNññîˆÇÃóvëfÇéÊÇËèoÇ∑
		int nowstack = stack[stack.size() - 1];
		stack.pop_back();
		//ÉXÉ^ÉbÉNÇ©ÇÁéÊÇËèoÇµÇΩóvëfÇ…ó◊ê⁄Ç∑ÇÈóvëfÇéÊìæ
		int neighbortri = _element[nowstack].neighbor[0];
		//ó◊ê⁄Ç∑ÇÈéOäpå`óvëfÇ™ë∂ç›Ç∑ÇÈÇ∆Ç´
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

				//ã´äEï”îªíË
				if (_element[nowstack].node[2] == _nodenumm1 || _element[nowstack].node[2] == _nodenump1) {
					_element[nowstack].side[1] = true;
					_element[neighbortri].side[2] = true;
				}

				_element[nowtri].getangle(_node);
				_element[neighbortri].getangle(_node);
				
				stack.push_back(nowstack);
				stack.push_back(neighbortri);
			}
		}
	}
}


//*****************************************************************************
//ã´äEÇÃê∂ê¨
//*****************************************************************************
void DelaunayClass::getboundary(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary) {
	getelement(_node, _element, -2, _boundary.nodelist[0], -2);
	for (int i = 1; i < _boundary.nodelist.size(); i++) {
		getelement(_node, _element, _boundary.nodelist[(i + 1) % _boundary.nodelist.size()], _boundary.nodelist[i], _boundary.nodelist[i - 1]);
	}
}


//*****************************************************************************
//ã´äEäOïîÇÃóvëfÇñ≥å¯âª
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
//óvëfóvëfä‘ó◊ê⁄ä÷åWÇÃèCïú
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
//ì‡ïîì_ÇÃê∂ê¨
//*****************************************************************************
void DelaunayClass::getinternalelement(vector<NodeClass> &_node, vector<ElementClass> &_element) {
	//ñ êœÇÃëÂÇ´Ç¢óvëfÇèdêSÇ≈ï™äÑ
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

		//í«â¡ÇµÇΩêﬂì_Ç≈óvëfï™äÑ
		getelement(_node, _element, -2, _node.size() - 1, -2);
	}
	//Ç–Ç∏Ç›ÇÃëÂÇ´Ç¢óvëfÇï™äÑ
	for (int i = 0; i < ADDITIONALNODENUM2; i++) {
		//í«â¡Ç∑ÇÈêﬂì_Çê∂ê¨
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

		//í«â¡ÇµÇΩêﬂì_Ç≈óvëfï™äÑ
		getelement(_node, _element, -2, _node.size() - 1, -2);
	}
}


//*****************************************************************************
//ñ≥å¯âªÇ≥ÇÍÇΩóvëfÇçÌèú
//*****************************************************************************
void DelaunayClass::deleteelement(vector<ElementClass> &_element) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		if (_element[i].active == false) {
			_element.erase(_element.begin() + i);
		}
	}
}