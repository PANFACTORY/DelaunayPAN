//*****************************************************************************
//Title		:Delaunay SolverClass
//Purpose	:Solver for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/27�`
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************


#include "DelaunayClass.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>


#define ADDITIONALNODENUM0	100000			//�ӂ̒����ɂ��ו�����


DelaunayClass::DelaunayClass() {}


DelaunayClass::~DelaunayClass() {}


//*****************************************************************************
//DelaunayTriangulation�̈�A�̏���
//*****************************************************************************
void DelaunayClass::delaunaymain(vector<NodeClass> &_node, vector<ElementClass> &_element, vector<BoundaryClass> &_boundary, double _maxsize, int _laplaciannum) {
	//----------SuperTriangle�̐���----------
	getsupertriangle(_node, _element);

	//----------���E�̐���----------
	for (int i = 0; i < _boundary.size(); i++) {
		getboundary(_node, _element, _boundary[i]);
	}
	for (int i = 0; i < _boundary.size(); i++) {
		deactivate(_node, _element, _boundary[i]);
	}

	//----------�s�v�ȗv�f�̍폜----------
	deletesupertriangle(_node, _element);
	deleteelement(_element);

	//----------�v�f�v�f�ԗאڊ֌W���C��----------
	sortelement(_element);

	if (_maxsize > 0) {
		//----------�ǉ��Őߓ_��z�u----------
		getinternalelement(_node, _element, _maxsize);

		//----------Laplacian�@�ɂ��ߓ_�̏C��----------
		laplacian(_node, _element, _laplaciannum);
	}
}


//*****************************************************************************
//�v�f���ɓ_�������Ƃ�
//*****************************************************************************
void DelaunayClass::getelementin(vector<NodeClass> &_node, vector<ElementClass> &_element, int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1) {
	vector<int> stack;
	vector<ElementClass> tmpelement(3);

	tmpelement[0].side[0] = _element[_nowtri].side[0];
	tmpelement[1].side[0] = _element[_nowtri].side[1];
	tmpelement[2].side[0] = _element[_nowtri].side[2];

	//���E�Ӕ���
	if (_element[_nowtri].node[0] == _nodenumm1 || _element[_nowtri].node[0] == _nodenump1) {
		tmpelement[1].side[1] = true;
		tmpelement[2].side[2] = true;
	}
	if (_element[_nowtri].node[1] == _nodenumm1 || _element[_nowtri].node[1] == _nodenump1) {
		tmpelement[2].side[1] = true;
		tmpelement[0].side[2] = true;
	}
	if (_element[_nowtri].node[2] == _nodenumm1 || _element[_nowtri].node[2] == _nodenump1) {
		tmpelement[0].side[1] = true;
		tmpelement[1].side[2] = true;
	}

	tmpelement[0].setnode(_nodenum, _element[_nowtri].node[1], _element[_nowtri].node[2]);
	tmpelement[0].setneighbor(_element[_nowtri].neighbor[0], _element.size(), _element.size() + 1);

	tmpelement[1].setnode(_nodenum, _element[_nowtri].node[2], _element[_nowtri].node[0]);
	tmpelement[1].setneighbor(_element[_nowtri].neighbor[1], _element.size() + 1, _nowtri);

	tmpelement[2].setnode(_nodenum, _element[_nowtri].node[0], _element[_nowtri].node[1]);
	tmpelement[2].setneighbor(_element[_nowtri].neighbor[2], _nowtri, _element.size());

	tmpelement[0].getangle(_node);
	tmpelement[1].getangle(_node);
	tmpelement[2].getangle(_node);

	for (int k = 0; k < 2; k++) {
		int neighbor = _element[_nowtri].neighbor[1 + k];
		if (neighbor >= 0) {
			_element[neighbor].neighbor[_element[neighbor].oppositenode(_nowtri)] = _element.size() + k;
		}
	}

	stack.push_back(_nowtri);
	stack.push_back(_element.size());
	stack.push_back(_element.size() + 1);

	_element[_nowtri].copy(tmpelement[0]);
	_element.push_back(tmpelement[1]);
	_element.push_back(tmpelement[2]);

	//.....�X���b�s���O.....
	swapping(_node, _element, stack, _nodenump1, _nodenumm1);
}


//*****************************************************************************
//�ӏ�ɓ_�������Ƃ�
//*****************************************************************************
void DelaunayClass::getelementon(vector<NodeClass> &_node, vector<ElementClass> &_element, int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1) {
	vector<int> stack;
	int nownode = _pos;
	int neitri = _element[_nowtri].neighbor[nownode];
	//�ӂ�����ŗאڗv�f�����݂���Ƃ�
	if (neitri != -1 && _element[neitri].active == true) {
		int neinode = _element[neitri].oppositenode(_nowtri);
		vector<ElementClass> tmptri(4);			//0:nowtri	1:neitri

		tmptri[0].setnode(_nodenum, _element[_nowtri].node[nownode], _element[_nowtri].node[(nownode + 1) % 3]);
		tmptri[0].setneighbor(_element[_nowtri].neighbor[(nownode + 2) % 3], neitri, _element.size());
		tmptri[0].setside(_element[_nowtri].side[(nownode + 2) % 3], false, false);

		tmptri[1].setnode(_nodenum, _element[_nowtri].node[(nownode + 1) % 3], _element[neitri].node[neinode]);
		tmptri[1].setneighbor(_element[neitri].neighbor[(neinode + 1) % 3], _element.size() + 1, _nowtri);
		tmptri[1].setside(_element[neitri].side[(neinode + 1) % 3], false, false);

		tmptri[2].setnode(_nodenum, _element[_nowtri].node[(nownode + 2) % 3], _element[_nowtri].node[nownode]);
		tmptri[2].setneighbor(_element[_nowtri].neighbor[(nownode + 1) % 3], _nowtri, _element.size() + 1);
		tmptri[2].side[0] = _element[_nowtri].side[(nownode + 1) % 3];

		tmptri[3].setnode(_nodenum, _element[neitri].node[neinode], _element[_nowtri].node[(nownode + 2) % 3]);
		tmptri[3].setneighbor(_element[neitri].neighbor[(neinode + 2) % 3], _element.size(), neitri);
		tmptri[3].side[0] = _element[neitri].side[(neinode + 2) % 3];

		int nei1 = _element[_nowtri].neighbor[(nownode + 1) % 3];
		if (nei1 != -1) {
			_element[nei1].neighbor[_element[nei1].oppositenode(_nowtri)] = _element.size();
		}

		int nei2 = _element[neitri].neighbor[(neinode + 2) % 3];
		if (nei2 != -1) {
			_element[nei2].neighbor[_element[nei2].oppositenode(neitri)] = _element.size() + 1;
		}

		if (tmptri[0].node[1] == _nodenumm1 || tmptri[0].node[1] == _nodenump1) {
			tmptri[2].side[1] = true;
			tmptri[0].side[2] = true;
		}
		if (tmptri[1].node[1] == _nodenumm1 || tmptri[1].node[1] == _nodenump1) {
			tmptri[0].side[1] = true;
			tmptri[1].side[2] = true;
		}
		if (tmptri[2].node[1] == _nodenumm1 || tmptri[2].node[1] == _nodenump1) {
			tmptri[2].side[2] = true;
			tmptri[3].side[1] = true;
		}
		if (tmptri[3].node[1] == _nodenumm1 || tmptri[3].node[1] == _nodenump1) {
			tmptri[3].side[2] = true;
			tmptri[1].side[1] = true;
		}

		if (_element[_nowtri].side[nownode] == true) {
			tmptri[0].side[1] = true;
			tmptri[1].side[2] = true;
			tmptri[2].side[2] = true;
			tmptri[3].side[1] = true;
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
		vector<ElementClass> tmptri(2);

		tmptri[0].setnode(_nodenum, _element[_nowtri].node[nownode], _element[_nowtri].node[(nownode + 1) % 3]);
		tmptri[0].setneighbor(_element[_nowtri].neighbor[(nownode + 2) % 3], -1, _element.size());
		tmptri[0].setside(_element[_nowtri].side[(nownode + 2) % 3], false, false);

		tmptri[1].setnode(_nodenum, _element[_nowtri].node[(nownode + 2) % 3], _element[_nowtri].node[nownode]);
		tmptri[1].setneighbor(_element[_nowtri].neighbor[(nownode + 1) % 3], _nowtri, -1);
		tmptri[1].setside(_element[_nowtri].side[(nownode + 1) % 3], false, false);

		int nei1 = _element[_nowtri].neighbor[(nownode + 1) % 3];
		if (nei1 != -1) {
			_element[nei1].neighbor[_element[nei1].oppositenode(_nowtri)] = _element.size();
		}

		if (tmptri[0].node[1] == _nodenumm1 || tmptri[0].node[1] == _nodenump1) {
			tmptri[1].side[1] = true;
			tmptri[0].side[2] = true;
		}
		if (tmptri[0].node[2] == _nodenumm1 || tmptri[0].node[2] == _nodenump1) {
			tmptri[0].side[1] = true;
		}
		if (tmptri[1].node[1] == _nodenumm1 || tmptri[1].node[1] == _nodenump1) {
			tmptri[1].side[2] = true;
		}

		if (_element[_nowtri].side[nownode] == true) {
			tmptri[0].side[1] = true;
			tmptri[1].side[2] = true;
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


//*****************************************************************************
//�X���b�s���O
//*****************************************************************************
void DelaunayClass::swapping(vector<NodeClass> &_node, vector<ElementClass> &_element, vector<int> &_stack, int _nodenump1, int _nodenumm1) {
	while (_stack.size() > 0) {
		//�X�^�b�N�����̗v�f�����o��
		int nowstack = _stack[_stack.size() - 1];
		_stack.pop_back();
		//�X�^�b�N������o�����v�f�ɗאڂ���v�f���擾
		int neighbortri = _element[nowstack].neighbor[0];
		//�אڂ���O�p�`�v�f�����݂���Ƃ�
		if (neighbortri >= 0 && _element[neighbortri].active == true) {
			int neighbornode = _element[neighbortri].oppositenode(nowstack);
			double r0 = _node[_element[nowstack].node[1]].distance(_node[_element[nowstack].node[2]]);
			double r1 = _node[_element[nowstack].node[0]].distance(_node[_element[neighbortri].node[neighbornode]]);
			if (r0 > r1
				&& _element[nowstack].inouton(_element[neighbortri].node[neighbornode], _node) == -1
				&& _element[neighbortri].inouton(_element[nowstack].node[0], _node) == -(neighbornode + 1)
				&& _element[nowstack].side[0] == false) {
				//cout << "!";

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

				//���E�Ӕ���
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
//SuperTriangle�̐���
//*****************************************************************************
void DelaunayClass::getsupertriangle(vector<NodeClass> &_node, vector<ElementClass> &_element) {
	//���_����ł����ꂽ�_��T��
	double rmax = 0.0;
	for (int i = 0; i < _node.size(); i++) {
		if (rmax < sqrt(pow(_node[i].x, 2.0) + pow(_node[i].y, 2.0))) {
			rmax = sqrt(pow(_node[i].x, 2.0) + pow(_node[i].y, 2.0));
		}
	}
	//rmax��1.5�{�̒��a�����~����ډ~�Ɏ��O�p�`�𐶐�
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
//SuperTriangle�̖�����
//*****************************************************************************
void DelaunayClass::deletesupertriangle(vector<NodeClass> &_node, vector<ElementClass> &_element) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				if (_element[i].node[j] == _node.size() - 1 - k) {
					_element[i].active = false;
					//_element.erase(_element.begin() + i);
					break;
				}
			}
		}
	}
}


//*****************************************************************************
//���E�̐���
//*****************************************************************************
void DelaunayClass::getboundary(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary) {
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
					if (_element[nowtri].neighbor[abs(pos) - 1] >= 0) {
						nowtri = _element[nowtri].neighbor[abs(pos) - 1];
					}
					else {
						cout << "Out of triangle Error!\n";
					}
				}
				//�v�f���ɂ��鎞
				else if (pos == 0) {
					//cout << "\nin->" << nowtri << "\t";
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
					//cout << "\non->" << nowtri << "\t";
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


//*****************************************************************************
//�v�f�𖳌���
//*****************************************************************************
void DelaunayClass::deactivate(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		int nodeorder[3];
		for (int j = 0; j < 3; j++) {
			nodeorder[j] = _boundary.order(_element[i].node[j]);
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


//*****************************************************************************
//�v�f�v�f�ԗאڊ֌W�̏C��
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
//�����_�̐���
//*****************************************************************************
void DelaunayClass::getinternalelement(vector<NodeClass> &_node, vector<ElementClass> &_element, double _maxside) {
	//�ӂ̒����v�f�𕪊�
	for (int i = 0; i < ADDITIONALNODENUM0; i++) {
		//�ǉ�����ߓ_�𐶐�
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

		if (maxside < _maxside) {
			//cout << "\ni=" << i << "\t" << "maxside=" << maxside << "\n";
			break;
		}

		NodeClass addnode;
		addnode.x = (_node[_element[maxelement].node[(maxnode + 1) % 3]].x + _node[_element[maxelement].node[(maxnode + 2) % 3]].x) / 2.0;
		addnode.y = (_node[_element[maxelement].node[(maxnode + 1) % 3]].y + _node[_element[maxelement].node[(maxnode + 2) % 3]].y) / 2.0;
		_node.push_back(addnode);

		//�ǉ������ߓ_�ŗv�f����
		//cout << "\non-->" << maxelement << "\t";
		getelementon(_node, _element, maxelement, maxnode, -2, _node.size() - 1, -2);
	}
}


//*****************************************************************************
//���������ꂽ�v�f���폜
//*****************************************************************************
void DelaunayClass::deleteelement(vector<ElementClass> &_element) {
	for (int i = _element.size() - 1; i >= 0; i--) {
		if (_element[i].active == false || _element[i].check == false) {
			_element.erase(_element.begin() + i);
		}
	}
}


//*****************************************************************************
//Laplacian�@
//*****************************************************************************
void DelaunayClass::laplacian(vector<NodeClass> &_node, vector<ElementClass> &_element, int _maxnum) {
	vector<int> logstack;
	int logstacknum = 100;			//�����ߓ_�΂�����C�����Ȃ��悤�ɒ��߂ɏC���������̂��X�g�b�N
	for (int i = 0; i < _maxnum; i++) {
		//���p�ő�v�f�Ƃ��̗v�f�ԍ���T��
		double angmax = 0.0;
		int elemax = -1, nodemax = -1;
		for (int j = 0; j < _element.size(); j++) {
			for (int k = 0; k < 3; k++) {
				if (angmax < fabs(60.0 - _element[j].angle[k]) && _node[_element[j].node[k]].type == false) {
					angmax = fabs(60.0 - _element[j].angle[k]);
					elemax = j;
					nodemax = _element[j].node[k];
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
		vector<int> stack;
		int nowtri = elemax;
		do{
			stack.push_back(nowtri);
			nowtri = _element[nowtri].neighbor[(_element[nowtri].nodeorder(nodemax) + 1) % 3];
		}while (nowtri != elemax);

		//�V�������W�����߂�
		double xdash = 0.0, ydash = 0.0;
		for (int j = 0; j < stack.size(); j++) {
			xdash += _node[_element[stack[j]].node[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].x 
				+ _node[_element[stack[j]].node[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].x;
			ydash += _node[_element[stack[j]].node[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].y
				+ _node[_element[stack[j]].node[(_element[stack[j]].nodeorder(nodemax) + 1) % 3]].y;
		}
		xdash /= 2.0*(double)stack.size();
		ydash /= 2.0*(double)stack.size();

		_node[nodemax].x = xdash;
		_node[nodemax].y = ydash;

		//cout << "\nLaplacian->" << nodemax << "\t" << xdash << "\t" << ydash;

		//���p�̍Čv�Z
		for (int j = 0; j < stack.size(); j++) {
			_element[stack[j]].getangle(_node);
		}
	}
}