//
// Created by Daniel Adamovic on 26/09/23.
//
#include "coordinates.h"

double eucl_dist(coord& c1, coord& c2) {
    return sqrt(pow(c1.get_x()-c2.get_x(),2) + pow(c1.get_y() - c2.get_y(), 2));
}