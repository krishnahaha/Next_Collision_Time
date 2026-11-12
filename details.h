/*Aim: To define a cell list for the particles.*/

#ifndef DETAILS_H_INCLUDED
#define DETAILS_H_INCLUDED
#include<vector>

template<typename real>
class details
{
public:
    details() {}
    ~details() {}

    class details *prev, *next; //Each pointer occupies 8 bytes.
    class details *left, *right, *parent;
    class details *p2;
    real eventTime, tlast;
    int pi, counter, counter2, queue;
    uint8_t type, boxestraveldx, boxestraveldy, boxestraveldz;

    //std::vector<details> cell; //Occupies 24 bytes.
};

#endif // DETAILS_H_INCLUDED