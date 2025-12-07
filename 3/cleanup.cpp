#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

void deleteAllFiles(const std::string& folderPath) {
    try {
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (fs::is_regular_file(entry.path())) {
                fs::remove(entry.path());
            }
        }
        std::cout << "All files deleted successfully in: " << folderPath << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    std::string folderPath = "output/random_walk/distance";

    if (fs::exists(folderPath) && fs::is_directory(folderPath)) {
        deleteAllFiles(folderPath);
    } else {
        std::cerr << "Invalid folder path!" << std::endl;
    }
    folderPath = "output/random_walk/path";

    if (fs::exists(folderPath) && fs::is_directory(folderPath)) {
        deleteAllFiles(folderPath);
    } else {
        std::cerr << "Invalid folder path!" << std::endl;
    }
    folderPath = "output/random_walk/steps";

    if (fs::exists(folderPath) && fs::is_directory(folderPath)) {
        deleteAllFiles(folderPath);
    } else {
        std::cerr << "Invalid folder path!" << std::endl;
    }

    return 0;
}