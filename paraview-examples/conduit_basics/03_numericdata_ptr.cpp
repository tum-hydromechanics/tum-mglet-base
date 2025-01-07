#include <iostream>
#include "conduit.hpp"


using namespace conduit;

int main(){
    int vsize = 5;
    std::vector<float64> vals(vsize,0.0);
    for(int i=0;i<vsize;i++)
    {
        vals[i] = 3.1415 * i;
    }

    Node n;
    n["v_owned"] = vals;
    n["v_external"].set_external(vals);

    n.info().print();

    n.print();

    vals[1] = -1 * vals[1];
    n.print();



}