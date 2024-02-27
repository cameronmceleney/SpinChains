//
// Created by Cameron McEleney on 27/02/2024.
//

#ifndef SPINCHAINS_HASHSPECIALISATION_H
#define SPINCHAINS_HASHSPECIALISATION_H

// C++ Standard Libraries
#include <functional>

// C++ Third Party Library

// C++ User Libraries (General)
#include "CommonStructures.h"

namespace std {
    template<>
    struct hash<CommonStructures::Point3D> {
        size_t operator()(const CommonStructures::Point3D& point) const noexcept {
            // Extract coordinates directly if CommonStructures::Point3D provides public access
            // Otherwise, use appropriate getter methods if coordinates are private
            const auto& [x, y, z] = point.coordinates; // Assuming structured binding is valid

            // Initial hash values for each component
            size_t h1 = std::hash<int>{}(x);
            size_t h2 = std::hash<int>{}(y);
            size_t h3 = std::hash<int>{}(z);

            // Seed value for combined hash, can be any constant, here using a prime number
            size_t seed = 0x9e3779b9;

            // Combine hashes using a method to mix bits more thoroughly
            // This method is inspired by the approach used in Boost and the magic number from the golden ratio
            seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);

            return seed;
        }
    };
}

#endif //SPINCHAINS_HASHSPECIALISATION_H
