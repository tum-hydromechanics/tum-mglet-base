#include <iostream>
#include "conduit.hpp"


using namespace conduit;

int main(){
    int64 vals[4] = {100,200,300,400};

    Node n;
    n.set(vals,4);

    int64 *my_vals = n.as_int64_ptr();

    for(index_t i=0; i < 4; i++)
    {
        std::cout << "my_vals[" << i << "] = " << my_vals[i] << std::endl;
    }


    //part 2
    Node n2;

    // set with integer c++11 initializer list
    n2.set({100,200,300});
    n2.print();

    // assign with integer c++11 initializer list
    n2 = {100,200,300};
    n2.print();

    // set with floating point c++11 initializer list
    n2.set({1.0,2.0,3.0});
    n2.print();

    // assign with floating point c++11 initializer list
    n2 = {1.0,2.0,3.0};
    n2.print();



}