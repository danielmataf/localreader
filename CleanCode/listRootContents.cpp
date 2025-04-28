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
    //define target is "LD2" 
    std::string target = "LD2"; // Change this to "CuSn" or "CxC" as needed
    std::string file1 = std::string("/home/matamoros/RGDv7") + target + "_testlow.root";
    std::string file2 = std::string("/home/matamoros/Aprfull") + target + "_simlow.root";
    std::string file3 = std::string("/home/matamoros/Aprsymm") + target + "_simlow.root";

    listRootContents(file1);
    listRootContents(file2);
    listRootContents(file3);

    return 0;
}
