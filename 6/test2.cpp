#include <windows.h>
#include <iostream>
#include <vector>
#include <bitset>

void PrintNUMATopology() {
    DWORD length = 0;
    GetLogicalProcessorInformationEx(RelationNumaNode, nullptr, &length);

    std::vector<BYTE> buffer(length);
    SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* info =
        reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data());

    if (!GetLogicalProcessorInformationEx(RelationNumaNode, info, &length)) {
        std::cerr << "Failed to get NUMA info." << std::endl;
        return;
    }

    BYTE* ptr = buffer.data();
    BYTE* end = ptr + length;

    std::cout << "NUMA Topology:\n";

    while (ptr < end) {
        auto* nodeInfo = reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(ptr);
        if (nodeInfo->Relationship == RelationNumaNode) {
            USHORT nodeNumber = nodeInfo->NumaNode.NodeNumber;
            GROUP_AFFINITY affinity = nodeInfo->NumaNode.GroupMask;

            std::bitset<64> mask(affinity.Mask);

            std::cout << "NUMA Node " << nodeNumber << " contains logical processors: ";
            for (DWORD i = 0; i < 64; ++i) {
                if (mask.test(i)) {
                    std::cout << i << " ";
                }
            }
            std::cout << "\n";
        }

        ptr += nodeInfo->Size;
    }
}

int main() {
    PrintNUMATopology();
    return 0;
}
