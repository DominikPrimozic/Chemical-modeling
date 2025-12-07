#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

int main() {
    for (const auto& entry : fs::directory_iterator("output/LJ/C/t02")) {
        if (entry.path().extension() == ".bin") {
            fs::remove(entry.path());
        }
    }
}