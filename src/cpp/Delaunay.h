//*****************************************************************************
//Title		:src/cpp/Delaunay.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "Node.h"
#include "Element.h"
#include "Boundary.h"


namespace DelaunayPAN{
	class DelaunayClass
	{
public:
		DelaunayClass();
		~DelaunayClass();

		void delaunaymain(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, std::vector<BoundaryClass> &_boundary, double _maxsize, int _laplaciannum);

		void getsupertriangle(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element);
		void deletesupertriangle(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element);

		void getboundary(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, BoundaryClass _boundary);
		void deactivate(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, BoundaryClass _boundary);
		void deleteelement(std::vector<ElementClass> &_element);
		void sortelement(std::vector<ElementClass> &_element);
		void getinternalelement(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, double _maxside);
		void swapping(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, std::vector<int> &_stack, int _nodenump1, int _nodenumm1);

		void getelementin(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1);
		void getelementon(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1);

		void laplacian(std::vector<NodeClass> &_node, std::vector<ElementClass> &_element, int _maxnum);
	};
}