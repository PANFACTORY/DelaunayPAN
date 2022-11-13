/**
 * @file boundary.h
 * @author PANFACTORY
 * @brief
 * @version 0.1
 * @date 2022-11-13
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include <vector>

namespace DelaunayPAN {
class Boundary {
   public:
    Boundary() {}
    ~Boundary() {}
    Boundary(std::vector<int>& _nodelists, bool _type) {
        this->nodelists = _nodelists;
        this->type = _type;
    }

    std::vector<int> nodelists;  // list of nodes on boundary
    bool type;                   // type of boundary

    int order(int _nodenum) {
        for (int i = 0; i < this->nodelists.size(); i++) {
            if (this->nodelists[i] == _nodenum) {
                return i;
            }
        }
        return -1;
    }  // get order id on boundary
};
}  // namespace DelaunayPAN
