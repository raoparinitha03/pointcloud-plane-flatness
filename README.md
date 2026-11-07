# planeperc: planar surface extraction & flatness analysis (from scratch)

A small C++17 perception library plus a ROS 2 node that finds the planes in a 3D point cloud
(floors, walls, ceilings) and measures how far each surface deviates from flat.
**No PCL, no Open3D, no Eigen.** Every algorithm is implemented in this repo:

| Stage | File | Algorithm |
|---|---|---|
| Downsample | `voxel_grid.hpp` | hash-grid centroid per voxel, O(n) |
| Neighbour search | `kdtree.hpp` | k-d tree, median split on widest axis, kNN with a bounded max-heap |
| Denoise | `outlier_removal.hpp` | statistical outlier removal (mean k-NN distance, mu + k*sigma) |
| Plane detection | `ransac_plane.hpp` | RANSAC with adaptive iteration count, LS refinement, sequential multi-plane extraction |
| Plane fitting | `plane_fit.hpp` | covariance matrix + Jacobi eigen-solver (smallest eigenvector = normal) |
| Flatness | `flatness.hpp` | in-plane (u,v) grid of cell-median deviations -> robust bump/dent map |
| Glue | `pipeline.hpp` | end-to-end pipeline with per-stage timing |

## Build & run (no ROS needed)

```bash
cmake -S . -B build && cmake --build build -j
ctest --test-dir build --output-on-failure     # unit tests
./build/synthetic_eval build                   # accuracy + parameter sweep, writes build/sweep.csv, build/flatness_wall.ply
```
(or without CMake: `g++ -std=c++17 -O2 -Iinclude apps/synthetic_eval.cpp -o synthetic_eval`)

## Run on ROS 2

```bash
colcon build --packages-select planeperc && source install/setup.bash
ros2 run planeperc plane_node --ros-args -r points:=/your/depth/points -p voxel_leaf:=0.03 -p ransac_threshold:=0.01
# RViz: add PointCloud2 on /flatness_cloud, colour = deviation from plane (blue = -, red = +)
```
`plane_node.cpp` was written without a ROS 2 install available to compile against. Build it first and fix any typos.
Works with a Gazebo depth camera or 3D LiDAR plugin (remap `points`).

## Results on synthetic data (known ground truth)

Scene: 5 x 4 x 2.6 m room, 3000 pts/m2, 4 mm Gaussian depth noise, 3 % random outliers,
and a 30 mm Gaussian bump (sigma 15 cm) on one wall.
Default run (voxel 30 mm, RANSAC threshold 10 mm), see `docs/eval_output.txt`:

* normal error: **0.001-0.006 deg**, offset error **< 0.2 mm** for all 3 planes
* planted 30 mm bump measured as **+31 mm** (cell-median peak); flat surfaces show a +-7 mm noise floor
* ~550 ms per frame on 134k raw points on one CPU core (SOR dominates, see below)

![flatness map](docs/flatness_map.png)

Parameter sweep (3 seeds each; `docs/sweep_results.csv`):

* RANSAC threshold below ~2.5x the sensor noise (here < 10 mm) fragments planes into extra detections;
  a threshold far above the noise (50 mm) still finds the planes but widens the flatness band and inflates the bump estimate.
* Voxel size trades speed for accuracy: 100 mm voxels are ~15x faster than none but average away part of the bump (25 mm vs 30 mm).

## Honest limitations

* Validated on **synthetic** data and unit tests; not yet on real sensor data. Real depth noise grows with range, and real scans have
  occlusions and curved/clutter surfaces.
* Statistical outlier removal is the bottleneck (k-d tree build + 16-NN per point, single thread). Next steps: skip SOR or lower k
  (RANSAC is already robust), parallelise, or use a voxel-hashed neighbour search.
* Sequential RANSAC can merge two coplanar but disconnected surfaces; no connected-component check yet.
* No ICP / Kalman / sensor fusion yet (planned: track planes across frames with an EKF, refine pose with point-to-plane ICP).
