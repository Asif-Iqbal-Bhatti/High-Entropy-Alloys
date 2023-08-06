#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm> // For std::transform

//================================================================
// AUTHOR :: Asif Iqbal
// USAGE  :: TO reorder POSCAR ATomic positions based on atom types
//           ./argv[0] POSCAR Li Cl P S ... NTypesAtom
//================================================================

struct PairHash {
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// Helper function to convert a string to lowercase
std::string toLowerCase(const std::string& str) {
    std::string lowerStr = str;
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
    return lowerStr;
}

std::ofstream outputFile("reorder_POSCAR");

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_filename> <type1> <type2> ... <typeN>\n";
        return 1;
    }

    std::ifstream inputFile(argv[1]);
    if (!inputFile) {
        std::cerr << "Error opening input file: " << argv[1] << '\n';
        return 1;
    }

    int Ntypes = 0;
    std::unordered_map<std::string, int> elementN;
    std::vector<std::string> element;
    std::unordered_map<std::pair<std::string, int>, std::string, PairHash> line;
    

    std::string lineStr;
    int lineCounter = 0;
    while (std::getline(inputFile, lineStr)) {
        ++lineCounter;

        if (lineCounter <= 5) {
            std::cout << lineStr << '\n';
            outputFile << lineStr << '\n';
        } 
        else if (lineCounter == 6) {
            std::istringstream iss(lineStr);
            std::string token;
            while (iss >> token) {
                element.push_back(token);
            }
            Ntypes = element.size();
        } 
        else if (lineCounter == 7) {
            std::istringstream iss(lineStr);
            for (int i = 0; i < Ntypes; ++i) {
                int count;
                iss >> count;
                elementN[element[i]] = count;
            }
        } 
        else if (toLowerCase(lineStr).find("direct") != std::string::npos) {
            for (int i = 0; i < Ntypes; ++i) {
                const std::string& atom = element[i];
                for (int j = 0; j < elementN[atom]; ++j) {
                    std::getline(inputFile, lineStr);
                    line[{atom, j + 1}] = lineStr;
                }
            }
        }
    }

    inputFile.close();

    for (int i = 2; i < argc; ++i) {
        std::cout << argv[i] << ' ';
    }
    std::cout << '\n';

    for (int i = 2; i < argc; ++i) {
        std::cout << elementN[argv[i]] << ' ';
    }
    std::cout << '\n';
    
    std::cout << "Direct\n";
    for (int i = 2; i < argc; ++i) {
        const std::string& atom = argv[i];
        for (int j = 1; j <= elementN[atom]; ++j) {
            std::cout << line[{atom, j}] << '\n';
        }
    }

    
    if (!outputFile) {
        std::cerr << "Error creating output file: reorder_POSCAR\n";
        return 1;
    }

    for (int i = 2; i < argc; ++i) {
        outputFile << argv[i] << ' ';
    }
    outputFile << '\n';

    for (int j = 2; j < argc; ++j) {
        outputFile << elementN[argv[j]] << ' ';
    }
    outputFile << '\n';
    
    outputFile << "Direct\n";
    for (int i = 2; i < argc; ++i) {
        const std::string& atom = argv[i];
        for (int j = 1; j <= elementN[atom]; ++j) {
            outputFile << line[{atom, j}] << '\n';
        }
    }

    outputFile.close();
    std::cerr << "File has been written to reorder_POSCAR\n";
    return 0;
}
