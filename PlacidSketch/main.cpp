#include "parm.h"
#include "stage1.h"
#include "stage2.h"
#include "stage3.h"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>

using namespace std;
class PacketProcessor {

    vector<Packet> packets;
    vector<string> csvFiles;

public:
    // Load all CSV files from a folder
    bool loadDataFromFolder(const string& folderPath) {
        try {
            // Get all CSV files from the folder
            for (const auto& entry : filesystem::directory_iterator(folderPath)) {
                if (entry.is_regular_file() && entry.path().extension() == ".csv") {
                    csvFiles.push_back(entry.path().string());
                }
            }

            sort(csvFiles.begin(), csvFiles.end());

            for (uint32_t windowNumber = 0; windowNumber < csvFiles.size(); ++windowNumber) {
                if (!loadSingleCSVFile(csvFiles[windowNumber], windowNumber)) {
                    cout << "Failed to load file: " << csvFiles[windowNumber] << endl;
                    return false;
                }
            }

            cout << "Total windows: " << csvFiles.size() << ", Total packets: " << packets.size() << endl;
            return true;

        } catch (const filesystem::filesystem_error& e) {
            cout << "Error accessing folder: " << e.what() << endl;
            return false;
        }
    }

    bool loadSingleCSVFile(const string& filename, uint32_t windowNumber) {
        ifstream dataFile(filename);
        string line;
        getline(dataFile, line); // Skip header line

        while (getline(dataFile, line)) {
            if (line.empty()) continue;

            size_t firstComma = line.find(',');
            string fingerprint, quintuple;

            if (firstComma != string::npos) {
                quintuple = line.substr(0, firstComma);
                size_t secondComma = line.find(',', firstComma + 1);
                if (secondComma != string::npos) {
                    size_t thirdComma = line.find(',', secondComma + 1);
                    fingerprint = (thirdComma != string::npos) 
                        ? line.substr(secondComma + 1, thirdComma - secondComma - 1)
                        : line.substr(secondComma + 1);
                } else {
                    fingerprint = line.substr(firstComma + 1);
                }
            } else {
                fingerprint = line;
            }

            if (!fingerprint.empty()) {
                packets.emplace_back(fingerprint.c_str(), quintuple.empty() ? nullptr : quintuple.c_str(), windowNumber);
            }
        }
        return true;
    }

    const vector<Packet>& getPackets() const { return packets; }
};


class PlacidSketch {
private:
    Stage3Merger stage3;
    Stage1Filter stage1;
    Stage2Monitor stage2;

    uint32_t currentWindow = 0;

public:
    explicit PlacidSketch(size_t stage1MemoryBytes = STAGE1_MEMORY_BYTES, 
                         size_t stage2MemoryBytes = STAGE2_MEMORY_BYTES)
        : stage1(stage1MemoryBytes), stage2(stage3, stage2MemoryBytes) {
    }

    void processPacket(const Packet& packet) {
        uint32_t windowSeq = packet.windowNumber;

        if (windowSeq != currentWindow) {
            stage1.resetBuckets(currentWindow);
            currentWindow = windowSeq;
        }

        if (stage1.processPacket(packet.flowID, windowSeq)) {
            stage2.processPotentialFlow(packet.flowID, windowSeq);
        }
    }
    void finalizeProcessing() {
        stage1.resetBuckets(currentWindow);
        stage3.finalize();
    }
};


int main() {
    cout << "PlacidSketch Stable Flow Detection" << endl;

    PacketProcessor dataLoader;
    string folderPath = "data"; // Change this to your data directory path
    
    if (!dataLoader.loadDataFromFolder(folderPath)) {
        return 1;
    }

    const auto& packets = dataLoader.getPackets();
    cout << "\n============== PlacidSketch Processing ==============" << endl;
    PlacidSketch sketch;
    for (const auto& packet : packets) {
        sketch.processPacket(packet);
    }
    sketch.finalizeProcessing();
    return 0;
}