# DH2323 Computer Graphics and Interaction
## Abstract
An accelerated raycast operation is proposed for the Unity game engine. High performance bounding volume hierarchies and efficient ray-box intersection tests are implemented in C++. A native plug-in with an easy-to-use API is constructed and integrated into Unity. Benchmarks are conducted that compare our implementation to Unity's built-in reference solution.

## References
> [1] Karras, T., 2012. Maximizing parallelism in the construction of BVHs, octrees, and k-d trees.
> In Proceedings of the Fourth ACM SIGGRAPH / Eurographics conference on High-Performance Graphics (EGGH-HPG’12).
> Eurographics Association, Goslar, DEU, 33–37.
> Available at: <https://devblogs.nvidia.com/wp-content/uploads/2012/11/karras2012hpg_paper.pdf>
   
> [2] Karras, T., 2012. Thinking Parallel, Part III: Tree Construction On The GPU | NVIDIA Developer Blog.
> NVIDIA Developer Blog.
> Available at: <https://devblogs.nvidia.com/thinking-parallel-part-iii-tree-construction-gpu/>

> [3] Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley. 2005. An efficient and robust ray-box intersection algorithm.
> In ACM SIGGRAPH 2005 Courses (SIGGRAPH ’05). Association for Computing Machinery, New York, NY, USA, 9–es.
> Available at: <https://doi.org/10.1145/1198555.1198748>, <https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans>

> [4] Santos, Artur L. dos, Alexandra Aníbal, Cidade Universitária, Veronica Teichrieb and Jorge Eduardo Falcao Lindoso. 2014.
> Review and Comparative Study of Ray Traversal Algorithms on a Modern GPU Architecture.
> Available at: <https://api.semanticscholar.org/CorpusID:44905165>