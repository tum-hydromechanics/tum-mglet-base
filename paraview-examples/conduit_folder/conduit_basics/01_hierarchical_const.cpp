#include <iostream>
#include "conduit.hpp"


using namespace conduit;

int main(){

    Node n;
    n["my"] = "data";
    n["a/b/c"] = "d";
    n["a"]["b"]["e"] = 64.0;
    n.print();

    std::cout << "total bytes: " << n.total_strided_bytes() << std::endl;
    std::cout<<"--"<<std::endl;

    /*
    Borrowing form JSON (and other similar notations), collections of named 
    nodes are called Objects and collections of unnamed nodes are called Lists,
    all other types are leaves that represent concrete data.
    */
    Node n1;
    n1["object_example/val1"] = "data";
    n1["object_example/val2"] = 10u;
    n1["object_example/val3"] = 3.1415;

    for(int i = 0; i < 5 ; i++ )
    {
        Node &list_entry = n1["list_example"].append();
        list_entry.set(i);
    }

    n1.print(); 

}