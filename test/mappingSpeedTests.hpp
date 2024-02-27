//
// Created by Cameron McEleney on 27/02/2024.
//

// C++ Standard Libraries
#include <iostream>
#include <map>
#include <unordered_map>
#include <random>
#include <vector>

// C++ User Libraries (General)
#include "../libs/CommonStructures.h"
#include "../libs/HashSpecialisation.h"

CommonStructures::Timer methodTimer;
bool randomInit;
bool showPrint;

void initialiseRandomSitesVector(std::vector<std::vector<std::vector<double>>>& grid, int& numSites,
                                 auto uni_rdm_real, auto rng) {
    int dimension = grid.size();
    for (int i = 0; i < numSites; ++i) {
        int x = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        int y = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        int z = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        grid[x][y][z] = uni_rdm_real(rng);
    }
}

void initialiseRandomSitesOM(std::map<CommonStructures::Point3D, double>& map, int numSites, int dimension,
                           std::uniform_real_distribution<double>& uni_rdm_real, std::mt19937& rng) {
    for (int i = 0; i < numSites; ++i) {
        int x = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        int y = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        int z = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        map[CommonStructures::Point3D(x, y, z)] = uni_rdm_real(rng);
    }
}

void initialiseRandomSitesUoM(std::unordered_map<CommonStructures::Point3D, double>& map, int numSites, int dimension,
                           std::uniform_real_distribution<double>& uni_rdm_real, std::mt19937& rng) {
    for (int i = 0; i < numSites; ++i) {
        int x = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        int y = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        int z = std::uniform_int_distribution<>(0, dimension - 1)(rng);
        map[CommonStructures::Point3D(x, y, z)] = uni_rdm_real(rng);
    }
}

void randomVector(int dimension_all, int num_computations, auto uni_rdm_real, auto uni_rdm_int, auto rng) {
    methodTimer.setName("Vector with Random Access");
    methodTimer.start(true);

    std::vector<std::vector<std::vector<double>>> nestedNestedVector(dimension_all, std::vector<std::vector<double>>(dimension_all, std::vector<double>(dimension_all)));

    // Populating the vector
    for (int i = 0; i < dimension_all; i++) {
        for (int j = 0; j < dimension_all; j++) {
            for (int k = 0; k < dimension_all; k++) {
                nestedNestedVector[i][j][k] = uni_rdm_real(rng); // Assigning a random value
            }
        }
    }

    for (int x = 0; x < dimension_all; ++x) {
        for (int y = 0; y < dimension_all; ++y) {
            for (int z = 0; z < dimension_all; ++z) {
                // Generate random coordinates within the dimensions and renormalise their sum
                double localSum = 0.0;
                for (int i = 0; i < num_computations; i++)
                    localSum += nestedNestedVector[uni_rdm_int(rng)][uni_rdm_int(rng)][uni_rdm_int(rng)];

                nestedNestedVector[x][y][z] = localSum / num_computations;
            }
        }
    }

    methodTimer.stop();
    methodTimer.print();

    for(auto& nestedVector : nestedNestedVector) {
        for(auto& v : nestedVector) {
            v.clear(); // Clears the innermost vector
            v.shrink_to_fit(); // Requests to reduce capacity of the innermost vector
        }
        nestedVector.clear(); // Clears the middle vector
        nestedVector.shrink_to_fit(); // Requests to reduce capacity of the middle vector
    }
    nestedNestedVector.clear(); // Clears the outer vector
    nestedNestedVector.shrink_to_fit(); // Requests to reduce capacity of the outer vector
}

void randomOrderedMap(int dimension_all, int num_computations, auto uni_rdm_real, auto uni_rdm_int, auto rng) {
    methodTimer.setName("Ordered Map with Random Access");
    methodTimer.start(true);

    std::map<CommonStructures::Point3D, double> exchangeContainerOrdered;

    for (int i = 0; i < dimension_all; i++) {
        for (int j = 0; j < dimension_all; ++j) {
            for (int k = 0; k < dimension_all; ++k) {
                exchangeContainerOrdered[CommonStructures::Point3D(i, j, k)] = uni_rdm_real(rng);;
            }
        }
    }

    if(exchangeContainerOrdered.size() >= 2) {
        for (auto& pair : exchangeContainerOrdered) {
            double localSum = 0.0;
            for (int i = 0; i < num_computations; i++) {
                // Generate random coordinates within the dimensions and renormalise their sum
                auto iter1 = exchangeContainerOrdered.find(CommonStructures::Point3D(uni_rdm_int(rng),
                                                                                     uni_rdm_int(rng),
                                                                                     uni_rdm_int(rng)));
                if (iter1 != exchangeContainerOrdered.end())
                    localSum += iter1->second;
                else
                    i--;
            }
            pair.second = localSum / num_computations;
        }
    }

    methodTimer.stop();
    methodTimer.print();

    exchangeContainerOrdered.clear();
    std::map<CommonStructures::Point3D, double>().swap(exchangeContainerOrdered);
}

void randomUnorderedMap(int dimension_all, int num_computations, auto uni_rdm_real, auto uni_rdm_int, auto rng) {
    methodTimer.setName("Unordered Map with Random Access");
    methodTimer.start(true);

    std::unordered_map<CommonStructures::Point3D, double> exchangeContainerUnordered;

    for (int i = 0; i < dimension_all; i++) {
        for (int j = 0; j < dimension_all; ++j) {
            for (int k = 0; k < dimension_all; ++k) {
                exchangeContainerUnordered[CommonStructures::Point3D(i, j, k)] = uni_rdm_real(rng);
            }
        }
    }

    if(exchangeContainerUnordered.size() >= 2) {
        for (auto& pair : exchangeContainerUnordered) {
            double localSum = 0.0;
            for (int i = 0; i < num_computations; i++) {
                // Generate random coordinates within the dimensions and renormalise their sum
                auto iter1 = exchangeContainerUnordered.find(CommonStructures::Point3D(uni_rdm_int(rng),
                                                                                     uni_rdm_int(rng),
                                                                                     uni_rdm_int(rng)));
                if (iter1 != exchangeContainerUnordered.end())
                    localSum += iter1->second;
                else
                    i--;
            }
            pair.second = localSum / num_computations;
        }
    }

    methodTimer.stop();
    methodTimer.print();

    exchangeContainerUnordered.clear();
    std::unordered_map<CommonStructures::Point3D, double>().swap(exchangeContainerUnordered);

}

void randomAccessTest(int container_size, int num_computations_input) {

    int dimension_all = container_size;
    int num_computations = num_computations_input;
    std::random_device rd;  // Initialises seed engine
    std::mt19937 rng(rd());  // Use Mersenne-Twister random-number engine

    double real_min = 0.0, real_max = 1.0;
    int int_min = 0.0, int_max = dimension_all - 1;
    std::uniform_real_distribution<double> uni_rdm_real(real_min, real_max);
    std::uniform_int_distribution<int> uni_rdm_int(int_min, int_max);

    randomVector(dimension_all, num_computations, uni_rdm_real, uni_rdm_int, rng);
    randomOrderedMap(dimension_all, num_computations, uni_rdm_real, uni_rdm_int, rng);
    randomUnorderedMap(dimension_all, num_computations, uni_rdm_real, uni_rdm_int, rng);

    std::exit(0);
}

void computeOnVector(int& dimension, int& numSites, auto uni_rdm_real, auto uni_rdm_int, auto rng) {



    std::vector<std::vector<std::vector<double>>> grid(dimension, std::vector<std::vector<double>>(dimension, std::vector<double>(dimension)));
    if (randomInit) {
        initialiseRandomSitesVector(grid, numSites, uni_rdm_real, rng);
    } else {
        // Populate the grid with initial values
        for (int x = 0; x < dimension; ++x) {
            for (int y = 0; y < dimension; ++y) {
                for (int z = 0; z < dimension; ++z) {
                    grid[x][y][z] = uni_rdm_real(rng);
                }
            }
        }
    }
    if (showPrint) {
        for ( int x = 0; x < dimension; ++x ) {
            for ( int y = 0; y < dimension; ++y ) {
                for ( int z = 0; z < dimension; ++z ) {
                    if ( grid[x][y][z] != 0.0 )
                        std::cout << "(" << x << ", " << y << ", " << z << "): " << grid[x][y][z] << " | ";
                }
            }
        }
    }


    methodTimer.setName("Vector with Neighbours Accessed");
    methodTimer.start(true);

    // Perform computation based on nearest neighbors
    for (int x = 0; x < dimension; ++x) {
        for (int y = 0; y < dimension; ++y) {
            for (int z = 0; z < dimension; ++z) {
                double sum = 0.0;
                int count = 0;

                for (int dx = -1; dx <= 1; dx += 2) {
                    int nx = x + dx;
                    if (nx >= 0 && nx < dimension) {
                        sum += grid[nx][y][z];
                        ++count;
                    }
                }

                for (int dy = -1; dy <= 1; dy += 2) {
                    int ny = y + dy;
                    if (ny >= 0 && ny < dimension) {
                        sum += grid[x][ny][z];
                        ++count;
                    }
                }

                for (int dz = -1; dz <= 1; dz += 2) {
                    int nz = z + dz;
                    if (nz >= 0 && nz < dimension) {
                        sum += grid[x][y][nz];
                        ++count;
                    }
                }

                // Update the cell based on the sum of its nearest neighbors
                grid[x][y][z] = sum / count;
            }
        }
    }

    methodTimer.stop();
    methodTimer.print();

    if (showPrint) {
        for ( int x = 0; x < dimension; ++x ) {
            for ( int y = 0; y < dimension; ++y ) {
                for ( int z = 0; z < dimension; ++z ) {
                    if ( grid[x][y][z] != 0.0 )
                        std::cout << "(" << x << ", " << y << ", " << z << "): " << grid[x][y][z] << " | ";
                }
            }
        }
    }
    std::cout << std::endl;

    for(auto& nestedVector : grid) {
        for(auto& v : nestedVector) {
            v.clear(); // Clears the innermost vector
            v.shrink_to_fit(); // Requests to reduce capacity of the innermost vector
        }
        nestedVector.clear(); // Clears the middle vector
        nestedVector.shrink_to_fit(); // Requests to reduce capacity of the middle vector
    }
    grid.clear(); // Clears the outer vector
    grid.shrink_to_fit(); // Requests to reduce capacity of the outer vector
}

void computeOnOrderedMap(int& dimension, int& numSites, auto uni_rdm_real, auto uni_rdm_int, auto rng) {

    std::map<CommonStructures::Point3D, double> map;

    if (randomInit) {
        initialiseRandomSitesOM(map, numSites, dimension, uni_rdm_real, rng);
    } else {
        // Populate the map with initial values
        for ( int x = 0; x < dimension; ++x ) {
            for ( int y = 0; y < dimension; ++y ) {
                for ( int z = 0; z < dimension; ++z ) {
                    map[CommonStructures::Point3D(x, y, z)] = uni_rdm_real(rng);
                }
            }
        }
    }

    if (showPrint) {
        for ( const auto &item: map ) {
            int x = item.first.x();
            int y = item.first.y();
            int z = item.first.z();
            std::cout << "(" << x << ", " << y << ", " << z << "): " << map[CommonStructures::Point3D(x, y, z)]
                      << " | ";
        }
    }

    // Temporary map to store updated values
    std::map<CommonStructures::Point3D, double> tempMap;

    methodTimer.setName("Ordered Map with Neighbours Accessed");
    methodTimer.start(true);

    // Perform computation based on nearest neighbors
    for (const auto& item : map) {
        double sum = 0.0;
        int count = 0;
        int x = item.first.x();
        int y = item.first.y();
        int z = item.first.z();

        for (int dx = -1; dx <= 1; dx += 2) {
            CommonStructures::Point3D neighbor(x + dx, y, z);
            if (map.find(neighbor) != map.end()) {
                sum += map[neighbor];
                ++count;
            }
        }

        for (int dy = -1; dy <= 1; dy += 2) {
            CommonStructures::Point3D neighbor(x, y + dy, z);
            if (map.find(neighbor) != map.end()) {
                sum += map[neighbor];
                ++count;
            }
        }

        for (int dz = -1; dz <= 1; dz += 2) {
            CommonStructures::Point3D neighbor(x, y, z + dz);
            if (map.find(neighbor) != map.end()) {
                sum += map[neighbor];
                ++count;
            }
        }

        // Update the cell based on the sum of its nearest neighbors
        tempMap[CommonStructures::Point3D(x, y, z)] = sum / std::max(1, count); // Avoid division by zero
    }

    // Copy the updated values from tempMap back to map
    map = std::move(tempMap);

    methodTimer.stop();
    methodTimer.print();

    if (showPrint) {
        for ( const auto &item: map ) {
            int x = item.first.x();
            int y = item.first.y();
            int z = item.first.z();
            std::cout << "(" << x << ", " << y << ", " << z << "): " << map[CommonStructures::Point3D(x, y, z)]
                      << " | ";
        }
    }

    std::cout << std::endl;


    map.clear();
    std::map<CommonStructures::Point3D, double>().swap(map);
}

void computeOnUnorderedMap(int& dimension, int& numSites, auto uni_rdm_real, auto uni_rdm_int, auto rng) {

    std::unordered_map<CommonStructures::Point3D, double> map;

    if (randomInit) {
        initialiseRandomSitesUoM(map, numSites, dimension, uni_rdm_real, rng);
    } else {
        // Populate the map with initial values
        for ( int x = 0; x < dimension; ++x ) {
            for ( int y = 0; y < dimension; ++y ) {
                for ( int z = 0; z < dimension; ++z ) {
                    map[CommonStructures::Point3D(x, y, z)] = uni_rdm_real(rng);
                }
            }
        }
    }

    if (showPrint) {
        for ( const auto &item: map ) {
            int x = item.first.x();
            int y = item.first.y();
            int z = item.first.z();
            std::cout << "(" << x << ", " << y << ", " << z << "): " << map[CommonStructures::Point3D(x, y, z)]
                      << " | ";
        }
    }

    // Temporary map to store updated values
    std::unordered_map<CommonStructures::Point3D, double> tempMap;

    methodTimer.setName("Unordered Map with Neighbours Accessed");
    methodTimer.start(true);

    // Perform computation based on nearest neighbors
    for (const auto& item : map) {
        double sum = 0.0;
        int count = 0;
        int x = item.first.x();
        int y = item.first.y();
        int z = item.first.z();

        for (int dx = -1; dx <= 1; dx += 2) {
            CommonStructures::Point3D neighbor(x + dx, y, z);
            if (map.find(neighbor) != map.end()) {
                sum += map[neighbor];
                ++count;
            }
        }

        for (int dy = -1; dy <= 1; dy += 2) {
            CommonStructures::Point3D neighbor(x, y + dy, z);
            if (map.find(neighbor) != map.end()) {
                sum += map[neighbor];
                ++count;
            }
        }

        for (int dz = -1; dz <= 1; dz += 2) {
            CommonStructures::Point3D neighbor(x, y, z + dz);
            if (map.find(neighbor) != map.end()) {
                sum += map[neighbor];
                ++count;
            }
        }

        // Update the cell based on the sum of its nearest neighbors
        tempMap[CommonStructures::Point3D(x, y, z)] = sum / std::max(1, count); // Avoid division by zero
    }

    // Copy the updated values from tempMap back to map
    map = std::move(tempMap);

    methodTimer.stop();
    methodTimer.print();

    if (showPrint) {
        for ( const auto &item: map ) {
            int x = item.first.x();
            int y = item.first.y();
            int z = item.first.z();
            std::cout << "(" << x << ", " << y << ", " << z << "): " << map[CommonStructures::Point3D(x, y, z)]
                      << " | ";
        }
    }

    std::cout << std::endl;


    map.clear();
    std::unordered_map<CommonStructures::Point3D, double>().swap(map);
}

void computeOnContainers(int dimension, int numSites, bool randomInitialisation = false, bool printOutput = false) {

    std::random_device rd;  // Initialises seed engine
    std::mt19937 rng(rd());  // Use Mersenne-Twister random-number engine

    double real_min = 0.0, real_max = 1.0;
    int int_min = 0, int_max = dimension - 1;
    std::uniform_real_distribution<double> uni_rdm_real(real_min, real_max);
    std::uniform_int_distribution<int> uni_rdm_int(int_min, int_max);

    randomInit = randomInitialisation;
    showPrint = printOutput;

    computeOnVector(dimension, numSites, uni_rdm_real, uni_rdm_int, rng);
    computeOnOrderedMap(dimension, numSites, uni_rdm_real, uni_rdm_int, rng);
    showPrint = true;
    computeOnUnorderedMap(dimension, numSites, uni_rdm_real, uni_rdm_int, rng);

    std::exit(0);

}

