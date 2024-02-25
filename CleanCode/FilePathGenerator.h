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

};

#endif // FILE_PATH_GENERATOR_H 