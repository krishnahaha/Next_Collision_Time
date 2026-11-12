/*Aim: To define a Neighbourlist for the particles*/

#ifndef NEIGHLIST_H_INCLUDED 
#define NEIGHLIST_H_INCLUDED

template<typename real>
class neighList{
public:
    std::vector<real> neigh;
    int nneigh;
};



#endif