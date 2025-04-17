#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>

void listRootContents(const std::string& filename) {
    TFile* file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return;
    }

    std::cout << "Contents of " << filename << ":\n";
    TIter next(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        std::string name = key->GetName();
        std::string cl   = key->GetClassName();
        std::cout << "  [" << cl << "] " << name << std::endl;
    }

    file->Close();
}

int main() {
    listRootContents("RinputFiles/Rhist_C2_RGD.root");
    return 0;
}
