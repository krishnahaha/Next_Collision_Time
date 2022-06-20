/*********************************************************************************************
 * Copyright (c) <2013>, <Santosh Ansumali@JNCASR>                                           *
 *  All rights reserved.                                                                     *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED       *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                  *
 *       Contributors:         Santosh Ansumali.                                             *
 *       Suggestions:          ansumali@jncasr.ac.in                                         *
 *       Bugs:                 ansumali@jncasr.ac.in                                         *
 *                                                                                           *
 *********************************************************************************************/

#ifndef _basic_matrix_ANSU_H_
#define  _basic_matrix_ANSU_H_

#include<iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include<iomanip>
using namespace std;

template<class dataType>
class matrix {
public:
  matrix(int nR, int nC, dataType scalar) {
    numRow = nR;
    numCol = nC;
    numElem = (nR)*(nC);
    data = new dataType[(nR)*(nC)];
    for(int row=0;row<nR;row++)
    for(int col=0;col<nC;col++)
    (*this)(row,col)= scalar;
}

~matrix()
  {
  if (this->numElem > 0)
  {
   for(int i =0;i<numRow;i++)
          data[i].~dataType();
        delete  data;
  }
  }

// copy constructor

matrix   (const matrix & tmp) {
  numRow  = tmp.numRow;
  numCol  = tmp.numCol;
  numElem = tmp.numElem;
  data = new dataType[numRow* numCol];
  std::memcpy(data , tmp.data, numRow* numCol* sizeof(dataType));
}

matrix & operator=(const matrix & tmp) {
    numRow  = tmp.numRow;
numCol  = tmp.numCol;
numElem = tmp.numElem;
if (&tmp == this){
          return *this;
      }

      for(int row=0;row<numRow;row++)
for(int col=0;col<numCol;col++)
(*this)(row,col)=tmp(row,col);


      return *this;
}
  dataType& operator()(int row, int col){
  return data[row * numCol + col];
}
const dataType operator()(int row, int col) const{

    return data[row * numCol + col];
}

 int nRow() const{return numRow;}

 int nCol()const{return numCol;}

  matrix    operator/(const  dataType  &second)
    {
       matrix    result( nRow(),  nCol(),0.0);
          for (int i=0; i < nRow(); i++)
            for (int  j=0; j <  nCol(); j++)
               result(i,j)= (*this)(i, j) / second;
        return result;
    }
matrix    operator*(const  dataType  &second)
    {
       matrix    result( nRow(),  nCol(),0.0);
          for (int i=0; i < nRow(); i++)
            for (int  j=0; j <  nCol(); j++)
               result(i,j)= (*this)(i, j) *second;
        return result;
    }
 template<class T>
friend std::ostream& operator << (std::ostream& os, const matrix<T> m);
template<class T>
friend matrix<T>  operator+(const matrix<T>  &first, const matrix<T> &second);
template<class T>
friend matrix<T>  operator*(const matrix<T>  &first, const matrix<T> &second);
dataType normL2(){
 dataType answer(0.0);
 for (int i=0; i <numRow; i++)
            for (int  j=0; j < numCol; j++)
                answer += (*this)(i,j)*(*this)(i,j);
return std::sqrt(answer);
}
private:
  dataType *data;
  int numElem;
  int numRow;
  int numCol;
};
   template<class dataType>
std::ostream& operator << (std::ostream& os, const matrix<dataType>  input) {
os << "(" << input.nRow() << "x" << input.nCol() << ") matrix:";
for (unsigned i= 0; i< input.nRow(); i++)
{ os << std::endl; // new row
for (unsigned j= 0; j< input.nCol(); j++)
os << std::setw(11) << std::right << input(i,j) << " " ;}
return os;}
template<class dataType>
  matrix<dataType>   operator+(const matrix<dataType>  &first, const matrix<dataType>  &second)
    {

        matrix<dataType>   result(first.nRow(), first.nCol(),0.0);
        for (int i=0; i <first.nRow(); i++)
            for (int  j=0; j < first.nCol(); j++)
                result(i, j) = first(i, j) + second(i, j);
        return result;
    }
template<class dataType>
  matrix<dataType>   operator-(const matrix<dataType>  &first, const matrix<dataType>  &second)
    {

        matrix<dataType>   result(first.nRow(), first.nCol(),0.0);
        for (int i=0; i <first.nRow(); i++)
            for (int  j=0; j < first.nCol(); j++)
                result(i, j) = first(i, j) - second(i, j);
        return result;
    }
template<class dataType>
  matrix<dataType>   operator*(const matrix<dataType>  &first, const matrix<dataType>  &second)
    {

        matrix<dataType>   result(first.nRow(), second.nCol(),0.0);
        for (int i=0; i <result.nRow(); i++)
            for (int  j=0; j < result.nCol(); j++){
result(i, j) = 0.0;
for (int  m=0; m < first.nCol(); m++)
                result(i, j) += first(i, m)*second(m, j);
}
        return result;
    }
template<class dataType>
  dataType   dotProduct(const matrix<dataType>  &first, const matrix<dataType>  &second){
dataType result = 0.0;
for (int i=0; i <first.nRow(); i++) {
result += first(i,0)*second(i,0);
  }
     return result;
  }
#endif // end of if
