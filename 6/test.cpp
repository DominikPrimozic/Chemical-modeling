#include <windows.h>
#include <iostream>
#include <vector>
#include <bitset>

void PrintProcessorInfo() {
    DWORD length = 0;
    GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &length);

    std::vector<BYTE> buffer(length);
    SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* info =
        reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data());

    if (!GetLogicalProcessorInformationEx(RelationProcessorCore, info, &length)) {
        std::cerr << "Failed to get processor info." << std::endl;
        return;
    }

    int coreIndex = 0;
    BYTE* ptr = buffer.data();
    BYTE* end = ptr + length;

    std::cout << "CPU Topology:\n";
    while (ptr < end) {
        auto* procInfo = reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(ptr);

        if (procInfo->Relationship == RelationProcessorCore) {
            GROUP_AFFINITY affinity = procInfo->Processor.GroupMask[0];
            std::bitset<64> mask(affinity.Mask);

            std::cout << "Physical Core " << coreIndex++ << " contains logical processors: ";
            for (DWORD i = 0; i < 64; ++i) {
                if (mask.test(i)) {
                    std::cout << i << " ";
                }
            }
            std::cout << "\n";
        }

        ptr += procInfo->Size;
    }
}

int main() {
    PrintProcessorInfo();
    return 0;
}