/**
 * @file node.h
 * @author PANFACTORY
 * @brief
 * @version 0.1
 * @date 2022-11-13
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

namespace DelaunayPAN {
template <class T>
class Node {
   public:
    Node() {
        this->x = T();
        this->y = T();
        this->isset = false;
        this->isonboundary = false;
    }
    ~Node() {}
    Node(T _x, T _y) {
        this->x = _x;
        this->y = _y;
        this->isset = false;
        this->isonboundary = false;
    }

    T x, y;             // cordinate vallue of Node
    bool isset;         // is set
    bool isonboundary;  // is on boundary

    T distance(const Node<T>& _node) const {
        return sqrt(pow(this->x - _node.x, 2.0) + pow(this->y - _node.y, 2.0));
    }
    T vecpro(const Node<T>& _node0, const Node<T>& _node1) const {
        return (_node0.x - this->x) * (_node1.y - this->y) -
               (_node0.y - this->y) * (_node1.x - this->x);
    }
    T innpro(const Node<T>& _node0, const Node<T>& _node1) const {
        return (_node0.x - this->x) * (_node1.x - this->x) +
               (_node1.y - this->y) * (_node0.y - this->y);
    }
};
}  // namespace DelaunayPAN
