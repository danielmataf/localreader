#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <chrono>
#include <thread>
#include "FilePathGenerator.h"
#include <stdexcept>

#include <unordered_set>
#include <algorithm>

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
    try {
        // Check if the directory exists and is not empty
        if (!fs::exists(directory)) {
            std::cerr << "Error: Directory does not exist: " << directory << std::endl;
            return;
        }

        filepaths.clear();  // Ensure filepaths vector is empty before populating

        // Iterate through directory and add all .hipo files to the vector
        int fileCount = 0;
        for (const auto& entry : fs::directory_iterator(directory)) {
            if (entry.path().extension() == ".hipo") {
                std::string filePath = entry.path().string();
                std::cout << "Found .hipo file: " << filePath << std::endl;  // Debugging line
                filepaths.push_back(filePath);
                fileCount++;
            }
        }

        // Post-iteration checks
        if (filepaths.empty()) {
            std::cerr << "No .hipo files found in the directory: " << directory << std::endl;
        } else {
            std::cout << "Total .hipo files added to vector: " << fileCount << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
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

// void FilePathGenerator::SnDir2Vector(const std::string& parentDirectory, std::vector<std::string>& filepaths) {
//        // this should go through all folders inside  '/volatile/clas12/dmat/gen/Sn/'
//        //std::cout << "pathfinding starts : " << parentDirectory << std::endl;
//        for (const auto& entry : fs::directory_iterator(parentDirectory)) {
//            if (entry.is_directory()) {
//                std::string subdirectory = entry.path().string();
//                std::cout << "subdirectory path: " << subdirectory << std::endl;
//                // check if the subdir has the file we want
//                std::string targetFile = subdirectory + "/sidis_mc-master/r_ttest3S.hipo";    //this is gonna be a problem... we need to change the name of the file in simu output 
//                //in the meantime another function LD2 specific has been created below
//                //wait, maybe not; To be continued  
//                if (fs::exists(targetFile)) {
//                    std::cout << "target file found: " << targetFile << std::endl;
//                    filepaths.push_back(targetFile);
//                }
//            }
//        }
//    }



void FilePathGenerator::SnDir2Vector(const std::string& parentDirectory, std::vector<std::string>& filepaths) {
    try {
        // Check if the parent directory exists before iterating
        if (!fs::exists(parentDirectory)) {
            throw std::runtime_error("Parent directory does not exist: " + parentDirectory);
        }

        for (const auto& entry : fs::directory_iterator(parentDirectory)) {
            if (entry.is_directory()) {
                std::string subdirectory = entry.path().string();
                std::cout << "Subdirectory path: " << subdirectory << std::endl;

                std::string targetFile = subdirectory + "/sidis_mc-master/r_ttest3S.hipo";
                
                if (fs::exists(targetFile)) {
                    std::cout << "Target file found: " << targetFile << std::endl;
                    filepaths.push_back(targetFile);
                } else {
                    std::cerr << "Target file not found in: " << subdirectory << std::endl;
                }
            }
        }

        // Extra validation to check if filepaths vector has elements before accessing it
        if (filepaths.empty()) {
            std::cerr << "No target files were found in the specified directory." << std::endl;
        } else {
            std::cout << "Total files found: " << filepaths.size() << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range error: " << e.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
    }
}
void FilePathGenerator::SnDir2VectorBis(const std::string& parentDirectory, std::vector<std::string>& filepaths) {
    try {
        // Ensure parent directory exists
        if (!fs::exists(parentDirectory)) {
            throw std::runtime_error("Parent directory does not exist: " + parentDirectory);
        }

        // Iterate over all subdirectories
        for (const auto& entry : fs::directory_iterator(parentDirectory)) {
            if (entry.is_directory()) {
                std::string subdirectory = entry.path().string();
                std::cout << "Subdirectory path: " << subdirectory << std::endl;

                std::string targetFile = subdirectory + "/sidis_mc-master/bis_ttest3S.hipo";
                
                // Check if the target file exists
                if (fs::exists(targetFile)) {
                    std::cout << "Target file found: " << targetFile << std::endl;
                    filepaths.push_back(targetFile);
                } else {
                    std::cerr << "Target file not found in: " << subdirectory << std::endl;
                }
            }
        }

        // Validation check: Ensure at least one file was found
        if (filepaths.empty()) {
            std::cerr << "No target files (bis_ttest3S.hipo) were found in the specified directory." << std::endl;
        } else {
            std::cout << "Total files found: " << filepaths.size() << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range error: " << e.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
    }
}


//
void FilePathGenerator::SnDir2VectorThird(const std::string& parentDirectory, std::vector<std::string>& filepaths) {
    try {
        // Ensure parent directory exists
        if (!fs::exists(parentDirectory)) {
            throw std::runtime_error("Parent directory does not exist: " + parentDirectory);
        }

        // Iterate over all subdirectories
        for (const auto& entry : fs::directory_iterator(parentDirectory)) {
            if (entry.is_directory()) {
                std::string subdirectory = entry.path().string();
                std::cout << "Subdirectory path: " << subdirectory << std::endl;

                std::string targetFile = subdirectory + "/sidis_mc-master/third_ttest3S.hipo";
                
                // Check if the target file exists
                if (fs::exists(targetFile)) {
                    std::cout << "Target file found: " << targetFile << std::endl;
                    filepaths.push_back(targetFile);
                } else {
                    std::cerr << "Target file not found in: " << subdirectory << std::endl;
                }
            }
        }

        // Validation check: Ensure at least one file was found
        if (filepaths.empty()) {
            std::cerr << "No target files (bis_ttest3S.hipo) were found in the specified directory." << std::endl;
        } else {
            std::cout << "Total files found: " << filepaths.size() << std::endl;
        }
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range error: " << e.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
    }
}



void FilePathGenerator::pass1search(const std::string& reconRoot,
                                    std::vector<std::string>& filepaths) {
    try {
        std::cout << "[pass1search] Start. Root: " << reconRoot << std::endl;

        // Step 1: reset output vector
        filepaths.clear();
        std::cout << "[step 1] Cleared output vector." << std::endl;

        // Step 2: validate directory
        if (!fs::exists(reconRoot) || !fs::is_directory(reconRoot)) {
            std::cerr << "[error] Not a directory: " << reconRoot << std::endl;
            return;
        }
        std::cout << "[step 2] Verified directory exists and is accessible." << std::endl;

        int fileCountRoot = 0;
        int subdirVisited = 0;
        int fileCountSub  = 0;

        // Step 3: pick up .hipo files directly under .../dst/recon/
        std::cout << "[step 3] Scanning root for .hipo files..." << std::endl;
        for (const auto& entry : fs::directory_iterator(reconRoot)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                const std::string filePath = entry.path().string();
                std::cout << "  [+] root .hipo: " << filePath << std::endl;
                filepaths.push_back(filePath);
                ++fileCountRoot;
            }
        }
        std::cout << "[step 3] Found " << fileCountRoot << " .hipo files at root." << std::endl;

        // Step 4: for each immediate subdirectory (e.g., 018530/), collect .hipo files
        std::cout << "[step 4] Scanning subdirectories for .hipo files..." << std::endl;
        for (const auto& sub : fs::directory_iterator(reconRoot)) {
            if (!sub.is_directory()) continue;

            ++subdirVisited;
            std::cout << "  [dir] Entering: " << sub.path().string() << std::endl;

            for (const auto& f : fs::directory_iterator(sub.path())) {
                if (f.is_regular_file() && f.path().extension() == ".hipo") {
                    const std::string filePath = f.path().string();
                    std::cout << "    [+] subdir .hipo: " << filePath << std::endl;
                    filepaths.push_back(filePath);
                    ++fileCountSub;
                }
            }
        }
        std::cout << "[step 4] Visited " << subdirVisited
                  << " subdirectories; found " << fileCountSub
                  << " .hipo files in subdirs." << std::endl;

        // Step 5: final report (same style as Files2Vector)
        const int total = fileCountRoot + fileCountSub;
        if (filepaths.empty()) {
            std::cerr << "[result] No .hipo files found under: " << reconRoot << std::endl;
        } else {
            std::cout << "[result] Total .hipo files added to vector: " << total << std::endl;
        }

        std::cout << "[pass1search] Done." << std::endl;
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "[filesystem error] " << e.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "[error] " << e.what() << std::endl;
    }
}


// Function to display progress with a percentage and a loading bar
//void FilePathGenerator::displayProgress(int current, int total) {
//    int barWidth = 50;
//    int progress = (current * barWidth) / total;
//    int percentage = (current * 100) / total;
//
//    std::cout << "\r[";
//    for (int i = 0; i < barWidth; ++i) {
//        if (i < progress) {
//            std::cout << "=";
//        } else if (i == progress) {
//            std::cout << ">";
//        } else {
//            std::cout << " ";
//        }
//    }
//    std::cout << "] " << percentage << "%";
//
//    std::cout.flush();
//}

void FilePathGenerator::displayProgress(int current, int total) {
    static int lastPercentage = -1; // Store the last displayed percentage
    int barWidth = 50;
    int percentage = (current * 100) / total;

    // Only update the display if the percentage has changed
    if (percentage > lastPercentage) {
        lastPercentage = percentage;

        int progress = (percentage * barWidth) / 100;

        std::cout << "\r[";
        for (int i = 0; i < barWidth; ++i) {
            if (i < progress) {
                std::cout << "=";
            } else if (i == progress) {
                std::cout << ">";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] " << percentage << "%";

        std::cout.flush();
    }
}


// New method to process events and show the progress
void FilePathGenerator::progressEvents(int totalEvents) {
    int counter = 0;

    // Example event processing loop
    for (int i = 0; i < totalEvents; ++i) {
        // Replace this with your actual event processing logic
        // std::optional<Event> event = ProcessNextEvent();
        
        // Simulate processing time
        std::this_thread::sleep_for(std::chrono::milliseconds(10)); // Remove or adjust as necessary

        // Update progress
        displayProgress(i + 1, totalEvents);
    }

    std::cout << "\nProcessing completed.\n";
}
