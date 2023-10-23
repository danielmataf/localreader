
int main() {
    std::vector<std::string> filenames;
    std::string line;

    // Open the text file containing file paths
    std::ifstream file("file_paths.txt");

    if (file.is_open()) {
        while (std::getline(file, line)) {
                    //looping in the file egetting each line that is a string 
                    
            filenames.push_back(line);
                //adding each line of txt file to the list filenames (vectr of strings)
        }

        file.close();
    } else {
        std::cerr << "Error opening file_paths.txt" << std::endl;
        return 1; 
    }


    std::cout << "Hello world\n";
    //...
    //...
    //...
        //
    return 0;
}