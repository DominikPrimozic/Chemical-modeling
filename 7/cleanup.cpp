#include <iostream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: cleanup <path_to_pat_directory>\n";
        return 1;
    }

    std::string pat = argv[1];

    if (!fs::exists(pat)) {
        std::cerr << "Directory does not exist: " << pat << std::endl;
        return 1;
    }

    try {
        // Delete pos/ and vel/ subdirectories
        if (fs::exists(pat + "/pos")) {
            fs::remove_all(pat + "/pos");
            std::cout << "Deleted: " << pat + "/pos" << std::endl;
        }
        if (fs::exists(pat + "/vel")) {
            fs::remove_all(pat + "/vel");
            std::cout << "Deleted: " << pat + "/vel" << std::endl;
        }

        // Delete metadata.txt and thermo.txt
        if (fs::exists(pat + "/metadata.txt")) {
            fs::remove(pat + "/metadata.txt");
            std::cout << "Deleted: " << pat + "/metadata.txt" << std::endl;
        }
        if (fs::exists(pat + "/thermo.txt")) {
            fs::remove(pat + "/thermo.txt");
            std::cout << "Deleted: " << pat + "/thermo.txt" << std::endl;
        }

        // Optionally delete pat directory if now empty
        if (fs::is_empty(pat)) {
            fs::remove(pat);
            std::cout << "Deleted empty directory: " << pat << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Cleanup complete.\n";
    return 0;
}
