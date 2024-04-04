#ifndef FILE_PATH_GENERATOR_H
#define FILE_PATH_GENERATOR_H

#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>

namespace fs = std::filesystem;

class FilePathGenerator {
public:
    // Default constructor
    FilePathGenerator();

    void readFilePathsFromTextFile(const std::string& filename, std::vector<std::string>& filepaths);
    void generateFilePathsFromDirectory(const std::string& directory, const std::string& filename);
    void Files2Vector(const std::string& directory, std::vector<std::string>& filepaths);
    void ParDir2Vector(const std::string& parentDirectory, std::vector<std::string>& filepaths); 
    //takes the path of parentDirectory, it will look the directories inside of it, enter in thel and append hipofiles to the vector filepaths
    void SnDir2Vector(const std::string& parentDirectory, std::vector<std::string>& filepaths);
    //this takes the specific path to recover simus, since we have called as numbers ans the hipo files are called the same way. see detail in .cpp


};

#endif // FILE_PATH_GENERATOR_H 