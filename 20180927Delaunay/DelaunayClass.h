//*****************************************************************************
//Title		:Delaunay SolverClass
//Purpose	:Solver for Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/27
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************

#include <vector>

#include "NodeClass.h"
#include "ElementClass.h"
#include "BoundaryClass.h"

using namespace std;

#pragma once
class DelaunayClass
{
public:
	DelaunayClass();
	~DelaunayClass();

	void getsupertriangle(vector<NodeClass> &_node, vector<ElementClass> &_element);
	void deletesupertriangle(vector<NodeClass> &_node, vector<ElementClass> &_element);
	void getelement(vector<NodeClass> &_node, vector<ElementClass> &_element, int _nodenump1, int _nodenum, int _nodenumm1);
	void getboundary(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary);
	void deleteelement(vector<NodeClass> &_node, vector<ElementClass> &_element, BoundaryClass _boundary);
	void sortelement(vector<ElementClass> &_element);
	void getinternalelement(vector<NodeClass> &_node, vector<ElementClass> &_element);
};

