//*****************************************************************************
//Title     :src/cpp/wrapper.cpp
//Author    :Tanabe Yuta
//Date      :2019/10/29
//Copyright :(C)2019 TanabeYuta
//*****************************************************************************


#define _USE_MATH_DEFINES


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tuple>
#include "Delaunay.h"


namespace py = pybind11;
using namespace DelaunayPAN;


template<class T>
std::tuple<std::vector<Node<T> >, std::vector<Element<T> > > delaunaymain2(std::vector<Node<T> >& _nodes, std::vector<Boundary>& _boundaries, T _maxsize, int _laplaciannum) {
    std::vector<Element<T> > elements;
    
    //----------Generate SuperTriangle----------
    getsupertriangle(_nodes, elements);

    //----------Generate Boundary----------
    for (auto& boundary : _boundaries) {
        getboundary(_nodes, elements, boundary);
    }
    for (auto& boundary : _boundaries) {
        deactivate(_nodes, elements, boundary);
    }

    //----------Delete needless Elements----------
    deletesupertriangle(_nodes, elements);
    deleteelement(elements);

    //----------Sort Elements----------
    sortelement(elements);

    if (_maxsize > 0) {
        //----------subdivide----------
        getinternalelement(_nodes, elements, _maxsize);

        //----------Laplacian smoothing----------
        laplacian(_nodes, elements, _laplaciannum);
    }
    
    return std::make_tuple(_nodes, elements);
}


PYBIND11_MODULE(DelaunayPAN, m) {
    //----------Node class----------
    py::class_<Node<double> >(m, "Node", "Node class of DelaunayPAN")
        .def(py::init<double, double>())
        .def("distance", &Node<double>::distance)
        .def_readwrite("x", &Node<double>::x)
        .def_readwrite("y", &Node<double>::y);
    
    //----------Boundary class----------
    py::class_<Boundary>(m, "Boundary", "Boundary class of DelaunayPAN")
        .def(py::init<std::vector<int>&, bool>())
        .def_readwrite("nodelists", &Boundary::nodelists)
        .def_readwrite("type", &Boundary::type);
    
    //----------Element class----------
    py::class_<Element<double> >(m, "Element", "Element class of DelaunayPAN")
        .def_readwrite("nodes", &Element<double>::nodes)
        .def_readwrite("neighbors", &Element<double>::neighbors)
        .def_readwrite("sides", &Element<double>::sides);

    //----------Delaunay----------
    m.def("delaunaymain", &delaunaymain2<double>);
}