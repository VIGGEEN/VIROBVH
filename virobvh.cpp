/**
 * Created by viktorm on 2020-05-21.
 *
 * [1] Karras, T., 2012. Maximizing parallelism in the construction of BVHs, octrees, and k-d trees.
 * In Proceedings of the Fourth ACM SIGGRAPH / Eurographics conference on High-Performance Graphics (EGGH-HPG’12).
 * Eurographics Association, Goslar, DEU, 33–37.
 * Available at: <https://devblogs.nvidia.com/wp-content/uploads/2012/11/karras2012hpg_paper.pdf>
 *
 * [2] Karras, T., 2012. Thinking Parallel, Part III: Tree Construction On The GPU | NVIDIA Developer Blog.
 * NVIDIA Developer Blog.
 * Available at: <https://devblogs.nvidia.com/thinking-parallel-part-iii-tree-construction-gpu/>
 *
 * [3] Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley. 2005. An efficient and robust ray-box intersection algorithm.
 * In ACM SIGGRAPH 2005 Courses (SIGGRAPH ’05). Association for Computing Machinery, New York, NY, USA, 9–es.
 * Available at: <https://doi.org/10.1145/1198555.1198748>, <https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans>
 *
 * [4] Santos, Artur L. dos, Alexandra Aníbal, Cidade Universitária, Veronica Teichrieb and Jorge Eduardo Falcao Lindoso.
 * 2014. Review and Comparative Study of Ray Traversal Algorithms on a Modern GPU Architecture.
 * Available at: <https://api.semanticscholar.org/CorpusID:44905165>
 */

#include "virobvh.h"
#include "ZippIterator.hpp"
#include <morton.h>
#include <execution>
#include <cmath>

// This makes the compiler happy, not used in practice
bool operator<(const Bounds &, const Bounds &) {
    return false;
}

bounding_volume_hierarchy::bounding_volume_hierarchy(Bounds *bounds, int32_t count) :
        internal_children((count - 1) * 2),
        internal_children_leaf((count - 1) * 2),
        internal_parents(count),
        leaf_parents(count),
        morton(count),
        internal_bounds(count),
        bounds(bounds, bounds + count) {

    // The idea is to assign a Morton code for each primitive ... [1]
    std::for_each(std::execution::par_unseq, morton.begin(), morton.end(), [&](uint64_t &element) {
        const size_t i = &element - &morton.front();
        const Vector3 center = (this->bounds.data() + i)->center();
        morton[i] = (uint64_t) libmorton::morton3D_64_encode(
                (uint_fast32_t) center.x,
                (uint_fast32_t) center.y,
                (uint_fast32_t) center.z
        );
    });

    // ... sort the Morton codes ... [1]
    auto zip = Zip(morton, this->bounds);
    std::sort(std::execution::par_unseq, zip.begin(), zip.end());

}

void bounding_volume_hierarchy::build_radix_tree() {

    // for each internal node with index i ∈ [0,n−2] in parallel [1]
    std::for_each(std::execution::par_unseq, morton.begin(), morton.end() - 1, [&](const uint64_t &element) {

        // Calculate current index
        const int idx = (&element - &morton[0]);

        // Determine direction of the range (+1 or -1) [1]
        const int d = (_delta(idx, idx + 1) - _delta(idx, idx - 1)) > 0 ? 1 : -1;

        // Compute upper bound for the length of the range [1]
        const auto deltaMin = _delta(idx, idx - d);
        int lMax = 2;
        while (_delta(idx, idx + lMax * d) > deltaMin) {
            lMax <<= 2;
        }

        // Find the other end using binary search [1]
        int l = 0;
        for (int t = lMax / 2; t; t /= 2) {
            if (_delta(idx, idx + (l + t) * d) > deltaMin) {
                l += t;
            }
        }

        // j ← i + l · d [1]
        const int j = idx + l * d;

        // Find the split position using binary search [1]
        const int deltaNode = _delta(idx, j); // the distance of the prefix of i

        int split = 0;
        int step = l;
        do {
            // exponential decrease [2]
            step = (step + 1) >> 1;
            if (_delta(idx, idx + (split + step) * d) > deltaNode) {
                split += step;
            }
        } while (step > 1);

        // Determine where to split the range. [2]
        const int finalSplit = idx + split * d + std::min<int>(d, 0);

        // Select childA. [2]
        internal_children[idx * 2] = finalSplit;
        internal_children_leaf[idx * 2] = std::min<int>(idx, j) == finalSplit;
        if (!internal_children_leaf[idx * 2]) {
            internal_parents[finalSplit] = idx;
        } else {
            leaf_parents[finalSplit] = idx;
        }

        // Select childB. [2]
        internal_children[idx * 2 + 1] = finalSplit + 1;
        internal_children_leaf[idx * 2 + 1] = std::max<int>(idx, j) == finalSplit + 1;
        if (!internal_children_leaf[idx * 2 + 1]) {
            internal_parents[finalSplit + 1] = idx;
        } else {
            leaf_parents[finalSplit + 1] = idx;
        }

    });
}

void bounding_volume_hierarchy::build_bounding_volume() {
    std::vector<std::atomic_int> visited(internal_parents.size());
    std::for_each(std::execution::par_unseq, leaf_parents.begin(), leaf_parents.end(),
                  [&](const int &point) {

                      // Calculate current index
                      size_t idx = &point - &leaf_parents.front();

                      // Each thread starts from one leaf node ... [1]
                      idx = leaf_parents[idx];

                      int vcount = 0;
                      while (true) {

                          // We track how many threads have visited each internal node using atomic counters ... [1]
                          vcount = visited[idx].fetch_add(1);

                          if (vcount == 0) {
                              // ... the first thread terminates immediately ... [1]
                              break;
                          } else {

                              // ... while  the  second  one  gets  to  process  the  node. [1]
                              const bool left_leaf = internal_children_leaf[idx * 2];
                              const int left_node = internal_children[idx * 2];
                              const bool right_leaf = internal_children_leaf[idx * 2 + 1];
                              const int right_node = internal_children[idx * 2 + 1];

                              //TODO: Should be able to vectorize this easily
                              internal_bounds[idx].minx = std::min(
                                      left_leaf ? bounds[left_node].minx : internal_bounds[left_node].minx,
                                      right_leaf ? bounds[right_node].minx : internal_bounds[right_node].minx
                              );
                              internal_bounds[idx].miny = std::min(
                                      left_leaf ? bounds[left_node].miny : internal_bounds[left_node].miny,
                                      right_leaf ? bounds[right_node].miny : internal_bounds[right_node].miny
                              );
                              internal_bounds[idx].minz = std::min(
                                      left_leaf ? bounds[left_node].minz : internal_bounds[left_node].minz,
                                      right_leaf ? bounds[right_node].minz : internal_bounds[right_node].minz
                              );
                              internal_bounds[idx].maxx = std::max(
                                      left_leaf ? bounds[left_node].maxx : internal_bounds[left_node].maxx,
                                      right_leaf ? bounds[right_node].maxx : internal_bounds[right_node].maxx
                              );
                              internal_bounds[idx].maxy = std::max(
                                      left_leaf ? bounds[left_node].maxy : internal_bounds[left_node].maxy,
                                      right_leaf ? bounds[right_node].maxy : internal_bounds[right_node].maxy
                              );
                              internal_bounds[idx].maxz = std::max(
                                      left_leaf ? bounds[left_node].maxz : internal_bounds[left_node].maxz,
                                      right_leaf ? bounds[right_node].maxz : internal_bounds[right_node].maxz
                              );
                          }

                          // If current is root, return
                          if (idx == 0) break;

                          // ... walks up the tree using parent pointers ... [1]
                          idx = internal_parents[idx];
                      }
                  });
}

// Traversing a BVH ... [4]
bool bounding_volume_hierarchy::intersection(const Ray &r) {
    return _intersection(r, 0);
}

// Traversing a BVH ... [4]
void bounding_volume_hierarchy::intersection_batch(const std::vector<Ray> &r, uint8_t* out, int count) {
    std::for_each(std::execution::par_unseq, r.begin(), r.end(),
                  [&](const Ray &ray) {
                      size_t idx = &ray - &r.front();
                      out[idx] = _intersection(ray, 0);
                  });
}

// Using δ(i, j) to denote the length of the longest common prefix between keys ki and kj [1]
int bounding_volume_hierarchy::_delta(const int i, const int j) {

    // Make sure we are not out of range
    if (j < 0 || j >= morton.size()) return -1;

    const uint64_t ki = morton[i];
    const uint64_t kj = morton[j];

    // The case of duplicate Morton codes has to be handled explicitly ... [1]
    // We accomplish this by augmenting each key with a bit representation of its index ... [1]
    // Identical Morton codes => split the range in the middle. [2]
    if (ki == kj) {
        return (i + j) >> 1;
    }

    // Calculate the number of highest bits that are the same
    // for all objects, using the count-leading-zeros intrinsic. [2]
    return __builtin_clzll(ki ^ kj);
}

// Traversing a BVH ... [4]
bool bounding_volume_hierarchy::_intersection(Ray r, int currentIndex) {

    int left = internal_children[currentIndex * 2];
    bool left_leaf = internal_children_leaf[currentIndex * 2];
    const Bounds &left_bounds = left_leaf ? bounds[left] : internal_bounds[left];
    const bool left_intersect = _raybox_intersection(left_bounds, r);
    if (left_intersect) {
        if (left_leaf) {
            return true;
        }
    }

    int right = internal_children[currentIndex * 2 + 1];
    bool right_leaf = internal_children_leaf[currentIndex * 2 + 1];
    const Bounds &right_bounds = right_leaf ? bounds[right] : internal_bounds[right];
    const bool right_intersect = _raybox_intersection(right_bounds, r);
    if (right_intersect) {
        if (right_leaf) {
            return true;
        }
    }

    return (left_intersect ? _intersection(r, left) : false) ||
           (right_intersect ? _intersection(r, right) : false);
}

// An efficient and robust ray-box intersection algorithm [3]
bool bounding_volume_hierarchy::_raybox_intersection(const Bounds &b, Ray r) {
    double tmin = -INFINITY, tmax = INFINITY;
    //TODO: Should be able to vectorize this easily
    {
        double t1 = (b.minx - r.origin.x) * (1.f / r.direction.x);
        double t2 = (b.maxx - r.origin.x) * (1.f / r.direction.x);
        tmin = std::max(tmin, std::min(t1, t2));
        tmax = std::min(tmax, std::max(t1, t2));
    }
    {
        double t1 = (b.miny - r.origin.y) * (1.f / r.direction.y);
        double t2 = (b.maxx - r.origin.y) * (1.f / r.direction.y);
        tmin = std::max(tmin, std::min(t1, t2));
        tmax = std::min(tmax, std::max(t1, t2));
    }
    {
        double t1 = (b.minz - r.origin.z) * (1.f / r.direction.z);
        double t2 = (b.maxz - r.origin.z) * (1.f / r.direction.z);
        tmin = std::max(tmin, std::min(t1, t2));
        tmax = std::min(tmax, std::max(t1, t2));
    }

    return tmax > std::max(tmin, 0.0);
}
