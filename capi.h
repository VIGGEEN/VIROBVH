//
// Created by viktorm on 2020-06-03.
//

#ifndef VIROBVH_CAPI_H
#define VIROBVH_CAPI_H

#include "virobvh.h"

extern "C" {

[[maybe_unused]] bounding_volume_hierarchy *create(float *bounds, int count);

[[maybe_unused]] bool intersection(bounding_volume_hierarchy * bvh, float * ray);

[[maybe_unused]] void intersection_batch(bounding_volume_hierarchy *bvh, float *ray, uint8_t * out, int count);

[[maybe_unused]] void destroy(bounding_volume_hierarchy * bvh);

}

#endif //VIROBVH_CAPI_H
