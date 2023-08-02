#pragma once

#include <absl/strings/str_cat.h>

#include <geode/basic/types.h>
#include <geode/basic/uuid.h>

#include <remfracres/CellId.h>
#include <geode/mesh/core/solid_mesh.h>
#include <cstddef>
#include <functional>
#include <string>


namespace remfracres {

        /*!
         * ID of a facet within a BRep.
         */
        struct FacetId {
            geode::uuid blockId;
            geode::index_t polyhedronIndex;
            geode::local_index_t facetIndex;

            /*!
             * Return the ID of the cell to which belongs the current facet.
             */
            CellId cellId() const {
                return CellId{ blockId, polyhedronIndex };
            }

            /*!
             * Return the PolyhedronFacet to be used to identify the facet within its parent mesh.
             */
            geode::PolyhedronFacet geodeFacet() const {
                return geode::PolyhedronFacet{ polyhedronIndex, facetIndex };
            }

            bool operator==(const FacetId &other) const {
                return blockId == other.blockId && polyhedronIndex == other.polyhedronIndex && facetIndex == other.facetIndex;
            }

            bool operator!=(const FacetId &other) const {
                return !operator==(other);
            }

            bool operator<(const FacetId &other) const {
                return blockId < other.blockId || (blockId == other.blockId && (polyhedronIndex < other.polyhedronIndex ||
                    (polyhedronIndex == other.polyhedronIndex && facetIndex < other.facetIndex)));
            }

            bool operator<=(const FacetId &other) const {
                return blockId < other.blockId || (blockId == other.blockId && (polyhedronIndex < other.polyhedronIndex ||
                    (polyhedronIndex == other.polyhedronIndex && facetIndex <= other.facetIndex)));
            }

            bool operator>(const FacetId &other) const {
                return !operator<=(other);
            }

            bool operator>=(const FacetId &other) const {
                return !operator>(other);
            }

            std::string string() const {
                return absl::StrCat("(", blockId.string(), ", ", polyhedronIndex, ", ", facetIndex, ")");
            }

        };

}


namespace std {

    template <>
    struct hash<remfracres::FacetId> {
        std::size_t operator()(const remfracres::FacetId &facetId) const {
            std::size_t polyhedronIndexHash = std::hash<geode::index_t>()(facetId.polyhedronIndex) << (sizeof(geode::local_index_t) * 8);
            std::size_t facetIndexHash = std::hash<geode::local_index_t>()(facetId.facetIndex);
            return std::hash<geode::uuid>()(facetId.blockId) ^ polyhedronIndexHash ^ facetIndexHash;
        }
    };

}
