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
	class Delaunay
	{
public:
		Delaunay();
		~Delaunay();

		void delaunaymain(std::vector<Node> &_node, std::vector<Element> &_element, std::vector<Boundary> &_boundary, double _maxsize, int _laplaciannum);

		void getsupertriangle(std::vector<Node> &_node, std::vector<Element> &_element);
		void deletesupertriangle(std::vector<Node> &_node, std::vector<Element> &_element);

		void getboundary(std::vector<Node> &_node, std::vector<Element> &_element, Boundary _boundary);
		void deactivate(std::vector<Node> &_node, std::vector<Element> &_element, Boundary _boundary);
		void deleteelement(std::vector<Element> &_element);
		void sortelement(std::vector<Element> &_element);
		void getinternalelement(std::vector<Node> &_node, std::vector<Element> &_element, double _maxside);
		void swapping(std::vector<Node> &_node, std::vector<Element> &_element, std::vector<int> &_stack, int _nodenump1, int _nodenumm1);

		void getelementin(std::vector<Node> &_node, std::vector<Element> &_element, int _nowtri, int _nodenump1, int _nodenum, int _nodenumm1);
		void getelementon(std::vector<Node> &_node, std::vector<Element> &_element, int _nowtri, int _pos, int _nodenump1, int _nodenum, int _nodenumm1);

		void laplacian(std::vector<Node> &_node, std::vector<Element> &_element, int _maxnum);
	};
}