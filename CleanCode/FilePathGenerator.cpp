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

void FilePathGenerator::ParDir2Vector(const std::string& parentDirectory, std::vector<std::string>& filepaths) {
    //This function takes a folder, and goes inside all the directories inside to recover every hipo file and put it in a vector
    std::vector<std::string> allFolders;

    for (const auto& entry : fs::directory_iterator(parentDirectory)) {
        if (entry.is_directory()) {
            allFolders.push_back(entry.path().string());
            //folderpaths.push_back(entry.path().string()); // Store the folder path
            //std::cout<<entry.path().string()<<std::endl;
            std::cout<<"Folder path: "<< allFolders.back()<<std::endl;
            Files2Vector(entry.path().string(), filepaths);
        }
    }
}

 void FilePathGenerator::SnDir2Vector(const std::string& parentDirectory, std::vector<std::string>& filepaths) {
        // this should go through all folders inside  '/volatile/clas12/dmat/gen/Sn/'
        for (const auto& entry : fs::directory_iterator(parentDirectory)) {
            if (entry.is_directory()) {
                std::string subdirectory = entry.path().string();
                std::cout << "subdirectory path: " << subdirectory << std::endl;
                // check if the subdir has the file we want
                std::string targetFile = subdirectory + "/sidis_mc-master/r_ttest3S.hipo";    //this is gonna be a problem... we need to change the name of the file in simu output 
                //in the meantime another function LD2 specific has been created below
                //wait, maybe not; To be continued  
                if (fs::exists(targetFile)) {
                    std::cout << "target file found: " << targetFile << std::endl;
                    filepaths.push_back(targetFile);
                }
            }
        }
    }

