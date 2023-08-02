#pragma once

#include <absl/strings/str_cat.h>

#include <geode/basic/types.h>
#include <geode/basic/uuid.h>

#include <cstddef>
#include <functional>
#include <string>


namespace remfracres {


        /*!
         * ID of a cell within a BRep.
         */
        struct CellId {
            geode::uuid blockId;
            geode::index_t polyhedronIndex;

            bool operator==(const CellId &other) const {
                return blockId == other.blockId && polyhedronIndex == other.polyhedronIndex;
            }

            bool operator!=(const CellId &other) const {
                return !operator==(other);
            }

            bool operator<(const CellId &other) const {
                return blockId < other.blockId || (blockId == other.blockId && polyhedronIndex < other.polyhedronIndex);
            }

            bool operator<=(const CellId &other) const {
                return blockId < other.blockId || (blockId == other.blockId && polyhedronIndex <= other.polyhedronIndex);
            }

            bool operator>(const CellId &other) const {
                return !operator<=(other);
            }

            bool operator>=(const CellId &other) const {
                return !operator>(other);
            }

            std::string string() const {
                return absl::StrCat("(", blockId.string(), ", ", polyhedronIndex, ")");
            }

        };

}


namespace std {

    template <>
    struct hash<remfracres::CellId> {
        std::size_t operator()(const remfracres::CellId &cellId) const {
            return std::hash<geode::uuid>()(cellId.blockId) ^ std::hash<geode::index_t>()(cellId.polyhedronIndex);
        }
    };

}
