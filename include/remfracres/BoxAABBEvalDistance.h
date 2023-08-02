#pragma once

#include <absl/strings/str_cat.h>

#include <geode/basic/types.h>
#include <geode/basic/uuid.h>
#include <geode/geometry/distance.h>
#include <geode/geometry/bounding_box.h>
#include <geode/geometry/point.h>
#include <absl/container/flat_hash_map.h>
#include <absl/container/flat_hash_set.h>
#include <cstddef>
#include <functional>
#include <string>


namespace remfracres {

class BoxAABBEvalDistance3D
{
public:

	BoxAABBEvalDistance3D()
	{
	};
	BoxAABBEvalDistance3D(absl::Span< const geode::BoundingBox3D > bounding_boxes ): bounding_boxes_( bounding_boxes )
    {
    };
    ~BoxAABBEvalDistance3D()
    {
    	bounding_boxes_.empty();
    };

    std::tuple< double, geode::Point3D > operator()(const geode::Point3D& query,geode::index_t current_element_box ) const
    {
        const auto box_center =( bounding_boxes_[current_element_box].min()+ bounding_boxes_[current_element_box].max() ) / 2.;
        return std::make_tuple(geode::point_point_distance( box_center, query ), box_center );
    };

    void set_bounding_boxes(absl::Span< const geode::BoundingBox3D > bounding_boxes){ bounding_boxes_ = bounding_boxes;};

private:
    absl::Span< const geode::BoundingBox3D > bounding_boxes_;
};

}


