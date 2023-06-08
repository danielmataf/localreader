#include <iostream>
#include "TFile.h"
#include "TH1F.h"

void simple_root_example() {
    TFile file("example.root");
    TH1F hist("hist", "Example Histogram", 100, -5, 5);
    for (int i = 0; i < 100000; i++) {
        hist.Fill(3);
    }
    std::cout << "File contains " << file.GetNkeys() << " keys" << std::endl;
    hist.Draw();
}
