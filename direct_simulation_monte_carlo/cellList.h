/*Aim: To define a cell list for the particles.*/

#ifndef CELLLIST_H_INCLUDED
#define CELLLIST_H_INCLUDED
#include<vector>

template<typename real>
class cellList
{
    public:
        cellList() {}

        ~cellList() {}

        std::vector<int> cell;
};

#endif // CELLLIST_H_INCLUDED