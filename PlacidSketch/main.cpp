#include "parm.h"
#include "stage1.h"
#include "stage2.h"
#include "stage3.h"
#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class PacketProcessor {
    std::vector<Packet> packets;

public:
    bool loadDataFromFolder(const std::string& folderPath) {
        try {
            packets.clear();

            std::vector<std::string> datFiles;
            for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
                if (entry.is_regular_file() && entry.path().extension() == ".dat") {
                    datFiles.push_back(entry.path().string());
                }
            }

            std::sort(datFiles.begin(), datFiles.end());

            for (uint32_t windowNumber = 0; windowNumber < datFiles.size(); ++windowNumber) {
                if (!loadSingleDatFile(datFiles[windowNumber], windowNumber)) {
                    std::cerr << "Failed to load file: " << datFiles[windowNumber] << std::endl;
                    return false;
                }
            }
            return true;
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error accessing folder: " << e.what() << std::endl;
            return false;
        }
    }

    bool loadSingleDatFile(const std::string& filename, uint32_t windowNumber) {
        std::ifstream dataFile(filename, std::ios::binary);
        if (!dataFile.is_open()) {
            std::cerr << "Cannot open file: " << filename << std::endl;
            return false;
        }

        uint32_t flow32 = 0;
        while (dataFile.read(reinterpret_cast<char*>(&flow32), sizeof(flow32))) {
            packets.emplace_back(flow32, windowNumber);
        }

        return true;
    }

    const std::vector<Packet>& getPackets() const { return packets; }
};

class PlacidSketch {
private:
    Stage3Merger stage3;
    Stage1Filter stage1;
    Stage2Monitor stage2;
    uint32_t currentWindow = 0;
    std::vector<uint32_t> detectedStableFlows;

public:
    explicit PlacidSketch(size_t stage1MemoryBytes = STAGE1_MEMORY_BYTES,
                         size_t stage2MemoryBytes = STAGE2_MEMORY_BYTES)
        : stage1(stage1MemoryBytes), stage2(stage3, stage2MemoryBytes) {
        stage3.setDetectedFlowsCallback(&detectedStableFlows);
    }

    void processPacket(const Packet& packet) {
        uint32_t windowSeq = packet.windowNumber;

        if (windowSeq != currentWindow) {
            stage1.resetBuckets(currentWindow);
            currentWindow = windowSeq;
        }

        if (stage1.processPacket(packet.flowID, windowSeq)) {
            stage2.processPacket(packet.flowID, windowSeq);
            stage2.processPotentialFlow(packet.flowID, windowSeq);
        }
    }

    void finalizeProcessing() {
        stage1.resetBuckets(currentWindow);
        stage3.finalize();

        std::sort(detectedStableFlows.begin(), detectedStableFlows.end());
        detectedStableFlows.erase(std::unique(detectedStableFlows.begin(), detectedStableFlows.end()),
                                 detectedStableFlows.end());
    }

    const std::vector<uint32_t>& getDetectedStableFlows() const {
        return detectedStableFlows;
    }
};

int main(int argc, char* argv[]) {
    std::string folderPath = (argc > 1) ? argv[1] : ".";
    PacketProcessor dataLoader;

    if (!dataLoader.loadDataFromFolder(folderPath)) {
        std::cerr << "Usage: " << argv[0] << " [data_folder_path]" << std::endl;
        return 1;
    }

    PlacidSketch sketch;
    for (const auto& packet : dataLoader.getPackets()) {
        sketch.processPacket(packet);
    }
    sketch.finalizeProcessing();

    std::cout << sketch.getDetectedStableFlows().size() << std::endl;
    return 0;
}
