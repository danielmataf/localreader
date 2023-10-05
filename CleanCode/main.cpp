#include <iostream>
#include <vector>
#include "reader.h"
#include "Event.h"
#include "EventReader.h"
//#include "EventReader.cpp"


//compile :
//g++ -o analyse main.cpp `root-config --cflags --libs` -I/home/matamoros/Desktop/localreader/hipo4/


int main() {
    std::vector<std::string> filenames =	{"../../../files2read/r_eD-01.hipo","../../../files2read/r_eSn-01.hipo","../../../files2read/r_eSn-01.hipo","../../../files2read/r_eSn-02.hipo"};
    std::cout<< "Hello world \n";
    EventReader MC(filenames);
    Event test;
    for (int i=0; i<10; i++){
           std::cout<< "Event number  " << i << std::endl;

            test = MC.ProcessEventsInFile();
            
            
            //test.CalcKinematics();    //WHERE SHOULD THIS BE CALLED 
            test.calcAll();
            test.Print();


    }    


    return 0;
}
