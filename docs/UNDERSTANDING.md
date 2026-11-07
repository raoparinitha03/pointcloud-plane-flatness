# Reverse-engineering guide: know every line well enough to rebuild it

Read each section, close the file, and re-derive it on paper. If you can, you can defend it in an interview.

## 0. The big picture
Raw cloud (noisy, huge, full of junk) -> make it small and clean -> find the dominant flat surfaces -> measure how flat they are.
Voxel grid = fewer points. SOR = remove junk. RANSAC = find a plane despite junk. LS fit = polish it. Flatness map = the actual measurement.

## 1. Voxel grid (`voxel_grid.hpp`)
* Idea: key = (floor(x/leaf), floor(y/leaf), floor(z/leaf)); sum points per key; output centroids.
* Why `floor`, not `(int)`: truncation maps -0.5 and +0.5 into the same cell -> bug near the origin. (Unit test covers it.)
* Why it helps: speed (RANSAC cost is proportional to n) and noise: averaging m points cuts random noise by ~sqrt(m).
* Cost of too-large leaf: real features smaller than the leaf get averaged away (see sweep: 100 mm voxels underestimate the 30 mm bump).
* Complexity O(n) average (hash map).

## 2. K-d tree (`kdtree.hpp`)
* Build: pick the axis of widest spread, put the median point in a node, recurse left (smaller) and right (larger). `nth_element` is O(n) per level -> O(n log n).
* Query: go down the near side first; the far side can only hold a closer point if the distance to the splitting plane < current k-th best distance. That one pruning test is the whole trick.
* Max-heap of size k keeps the current k best; the top is the worst, so it is the pruning radius.
* Test: compare against brute force (`unit_tests.cpp::testKdTree`). Always test a clever algorithm against a dumb one.

## 3. Statistical outlier removal (`outlier_removal.hpp`)
* For each point, mean distance to its k nearest neighbours. Points on a surface have small values; flying pixels have big ones.
* Threshold = mean + stdMult * std of those means. Assumes the distribution is roughly Gaussian; heavy outlier contamination inflates sigma and weakens the filter (you saw leftover outliers in the first run).
* Query k+1 neighbours and skip the first: the nearest neighbour of a point is itself.

## 4. RANSAC plane (`ransac_plane.hpp`)
* Sample 3 points -> normal = (b-a) x (c-a), normalised; d = -n.a. Collinear points give a zero cross product -> skip.
* Score = number of points with |n.p + d| < threshold. Keep the best.
* **Adaptive iterations.** Probability a random 3-sample is all inliers = w^3 (w = inlier ratio). To be p-confident of at least one clean sample: N = log(1-p) / log(1-w^3). Example: w = 0.3, p = 0.99 -> N ~ 169. w = 0.5 -> N ~ 35. That is why the logs show 5-65 iterations, not 1000.
* **Threshold choice:** ~2.5-3 sigma of sensor noise. Too small -> real surface points rejected, planes fragment. Too big -> neighbouring surfaces merge and the bump disappears into the inliers.
* **Refinement:** the best 3-point plane is noisy. Fit all inliers with least squares, re-collect inliers, repeat (3 times). Accuracy goes from ~0.5 deg to ~0.002 deg.
* **Sequential extraction:** remove inliers of the best plane and run again, giving floor, wall, wall.
* Weakness: ignores connectivity (two separate coplanar patches become one plane); needs a threshold you must tune.

## 5. Least-squares plane via eigenvectors (`plane_fit.hpp`)
* Minimise sum of squared PERPENDICULAR distances. The plane passes through the centroid; the normal is the direction along which the points vary least.
* Covariance C = sum (p - c)(p - c)^T (3x3, symmetric). Its eigenvectors are the principal axes; the one with the smallest eigenvalue is the normal.
* **Jacobi method:** repeatedly apply a 2D rotation that zeroes one off-diagonal entry; after a few sweeps the matrix is diagonal, eigenvalues on the diagonal, accumulated rotations = eigenvectors. Angle: tan(2*theta) = 2a_pq / (a_pp - a_qq); the code uses the numerically stable `t` form from Numerical Recipes.
* Why not z = ax + by + c least squares? It minimises vertical error and fails for vertical planes (walls!). The eigenvector method has no preferred axis.
* Sign ambiguity: n and -n describe the same plane. `orientToward(viewpoint)` flips the normal to face the sensor so deviations have a consistent sign.

## 6. Flatness map (`flatness.hpp`)
* Problem: a bump taller than the RANSAC threshold is not an inlier, so inliers alone hide it. Fix: gather a wider band (+-5 cm) of points around the plane, restricted to the plane's own extent and excluding points owned by other planes.
* Problem 2: max |deviation| is useless; one outlier ruins it (first run: 40-50 mm "bumps" on flat walls). Fix: bin points into 5 cm cells in the plane's (u,v) coordinates, take the **median** per cell, report extremes of the medians. Robust statistics: a median ignores a few outliers.
* Noise floor: extremes over hundreds of noisy cells are ~+-7 mm even on a perfect plane, so anything much larger than that (30 mm) is real. Always know your noise floor before claiming a detection.

## 7. Questions an interviewer may ask (answer without looking)
1. Why RANSAC and not just least squares on everything? (LS is dragged by outliers: it has a 0 % breakdown point.)
2. How many iterations do you need and why? (formula above)
3. How did you choose the threshold? (sensor noise, sweep results)
4. Why a k-d tree? What is its worst case? (high dimensions / degenerate data -> O(n); voxel hashing is an alternative)
5. Why is your pipeline slow, and how would you speed it up? (SOR; profile numbers are printed per stage)
6. How would you add ICP / an EKF? (point-to-plane ICP between frames using the plane normals; EKF state = plane (n, d), predicted from IMU/odometry, updated by the RANSAC result)
7. What would break on a real sensor? (range-dependent noise, multipath, occlusion, non-planar clutter, calibration error)

## 8. Drills (do these; each takes 10-20 min)
1. Set RANSAC `maxIters` fixed at 1000 and disable the adaptive update. Compare iterations and time.
2. Set the threshold to 3 mm, then 100 mm; explain every change in the output.
3. Delete the `floor` -> `(int)` bug fix on purpose (use `static_cast<int>(p.x*inv)`) and watch the voxel unit test fail.
4. Replace the median with the mean in the flatness map and rerun with the outliers on.
5. Add 20 % outliers to `makeRoom`; at what point does SOR stop helping?
6. Implement point-to-plane ICP (one iteration = solve a 6x6 linear system) and test it on two rotated copies of the room.
