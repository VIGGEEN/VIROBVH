//
// Created by viktorm on 2020-06-03.
//
#include "capi.h"

[[maybe_unused]] bounding_volume_hierarchy *create(float *bounds, int count) {
    static_assert(sizeof(Bounds) == 6 * 4);
    count /= 6;
    bounding_volume_hierarchy* bvh = new bounding_volume_hierarchy(reinterpret_cast<Bounds *>(bounds), count);
    bvh->build_radix_tree();
    bvh->build_bounding_volume();
    return bvh;
}

[[maybe_unused]] bool intersection(bounding_volume_hierarchy *bvh, float *ray) {
    static_assert(sizeof(Ray) == 6 * 4);
    return bvh->intersection(*reinterpret_cast<Ray *>(ray));
}

[[maybe_unused]] void intersection_batch(bounding_volume_hierarchy *bvh, float *ray, uint8_t * out, int count) {
    static_assert(sizeof(Ray) == 6 * 4);
    count /= 6;
    std::vector<Ray> rays(reinterpret_cast<Ray *>(ray), reinterpret_cast<Ray *>(ray) + count);
    return bvh->intersection_batch(rays, out, count);
}

[[maybe_unused]] void destroy(bounding_volume_hierarchy *bvh) {
    delete bvh;
}
