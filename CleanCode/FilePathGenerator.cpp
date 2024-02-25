#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include "FilePathGenerator.h"

// Default constructor definition
FilePathGenerator::FilePathGenerator() {
}

void FilePathGenerator::generateFilePathsFromDirectory(const std::string& directory, const std::string& filename) {
    // Generate file paths automatically from a directory and write them to a text file
    std::ofstream file(filename, std::ios::trunc); // Open the file in truncate mode to overwrite it
    if (!file.is_open()) {
        std::cerr << "Unable to create file: " << filename << std::endl;
        return;
    }

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".hipo") {
            //file << entry.path() << '\n';
            //file << '"' << entry.path() << '"' << ',' << '\n'; // Add quotation marks and comma
        }
    }

    file.close();
}

void FilePathGenerator::readFilePathsFromTextFile(const std::string& filename, std::vector<std::string>& filepaths) {
    // Read file paths from a text file
    std::cout << "Current working directory: " << fs::current_path() << std::endl;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        filepaths.push_back(line);
    }

    file.close();
}

void FilePathGenerator::Files2Vector(const std::string& directory, std::vector<std::string>& filepaths) {
    // Generate file paths automatically from a directory and append them to the vector
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".hipo") {
            //std::cout << entry.path() << std::endl;
            filepaths.push_back(entry.path().string());
        }
    }
}