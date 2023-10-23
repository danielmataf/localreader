#include <iostream>
#include <fstream>
#include <filesystem>




// creating a function that gets a list of directories (= path to directories, here CuSn, LD2 and C)
// and gets also the name of an file (.txt ) in order to write in that text file
// th e function wshould get the content of all directories (file names) and write thel down in the txt output file 
// can be called anywhere ? 

void createFilePathsFile(const std::vector<std::string>& directories, const std::string& outputFilePath) {
    std::ofstream outputFile(outputFilePath);

    if (!outputFile.is_open()) {
        std::cerr << "unable to open  output file." << std::endl;
        return;
    }

    for (const std::string& directory : directories) {
            //creating an iterator in the directory list (vector of strings)
        std::filesystem::path path(directory);
            //definign a path from std lib
        if (std::filesystem::is_directory(path)) {
            for (const auto entry : std::filesystem::directory_iterator(path)) {
                if (entry.path().extension() == ".hipo") {
                    outputFile << entry.path() << "\n";
                }
            }
        }
    }

    outputFile.close();
}

int main() {
    std::vector<std::string> directories = {
        "/home/matamoros/Desktop/LumiScanDta/CuSn",
        "/home/matamoros/Desktop/LumiScanDta/LD2",
        "/home/matamoros/Desktop/LumiScanDta/C",
        // other direct can be added
    };

    std::string outputFilePath = "file_paths.txt";

    createFilePathsFile(directories, outputFilePath);

    std::cout << "File paths have been written to " << outputFilePath << std::endl;

    return 0;
}
