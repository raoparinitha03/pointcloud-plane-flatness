// ROS 2 wrapper around the from-scratch pipeline.
//   subscribes : points           (sensor_msgs/PointCloud2)  <- remap to your depth camera / 3D LiDAR topic
//   publishes  : flatness_cloud   (sensor_msgs/PointCloud2)  XYZRGB, colour = signed deviation from fitted plane
// Parameters  : voxel_leaf, ransac_threshold, min_inliers, max_planes, sor_k, sor_std, flat_band, flat_cell, deviation_range
//
// NOTE: this file was written without a ROS 2 install available to compile against.
//       Build it on your machine (ROS 2 Humble/Jazzy):  colcon build --packages-select planeperc
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>

#include "planeperc/pipeline.hpp"
#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"

class PlaneNode : public rclcpp::Node {
 public:
  PlaneNode() : Node("plane_extraction_node") {
    declare_parameter<double>("voxel_leaf", 0.03);
    declare_parameter<double>("ransac_threshold", 0.01);
    declare_parameter<int>("max_planes", 6);
    declare_parameter<int>("min_inliers", 200);
    declare_parameter<int>("sor_k", 16);
    declare_parameter<double>("sor_std", 1.5);
    declare_parameter<double>("flat_band", 0.05);
    declare_parameter<double>("flat_cell", 0.05);
    declare_parameter<double>("deviation_range", 0.03);

    sub_ = create_subscription<sensor_msgs::msg::PointCloud2>(
        "points", rclcpp::SensorDataQoS(), std::bind(&PlaneNode::onCloud, this, std::placeholders::_1));
    pub_ = create_publisher<sensor_msgs::msg::PointCloud2>("flatness_cloud", 1);
  }

 private:
  // Read x,y,z (FLOAT32) out of a PointCloud2 by hand - no pcl_conversions needed.
  static bool parseXYZ(const sensor_msgs::msg::PointCloud2& m, planeperc::Cloud& out) {
    int ox = -1, oy = -1, oz = -1;
    for (const auto& f : m.fields) {
      if (f.datatype != sensor_msgs::msg::PointField::FLOAT32) continue;
      if (f.name == "x") ox = static_cast<int>(f.offset);
      if (f.name == "y") oy = static_cast<int>(f.offset);
      if (f.name == "z") oz = static_cast<int>(f.offset);
    }
    if (ox < 0 || oy < 0 || oz < 0) return false;
    const size_t n = static_cast<size_t>(m.width) * m.height;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      const uint8_t* base = m.data.data() + i * m.point_step;
      float x, y, z;
      std::memcpy(&x, base + ox, 4);
      std::memcpy(&y, base + oy, 4);
      std::memcpy(&z, base + oz, 4);
      if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z)) out.push_back({x, y, z});
    }
    return true;
  }

  void onCloud(const sensor_msgs::msg::PointCloud2::ConstSharedPtr msg) {
    planeperc::Cloud raw;
    if (!parseXYZ(*msg, raw)) {
      RCLCPP_WARN_THROTTLE(get_logger(), *get_clock(), 2000, "PointCloud2 has no float32 x/y/z fields");
      return;
    }
    planeperc::PipelineParams prm;
    prm.voxelLeaf = get_parameter("voxel_leaf").as_double();
    prm.ransac.threshold = get_parameter("ransac_threshold").as_double();
    prm.maxPlanes = static_cast<int>(get_parameter("max_planes").as_int());
    prm.ransac.minInliers = static_cast<size_t>(get_parameter("min_inliers").as_int());
    prm.sorK = static_cast<size_t>(get_parameter("sor_k").as_int());
    prm.sorStd = get_parameter("sor_std").as_double();
    prm.flatBand = get_parameter("flat_band").as_double();
    prm.flatCell = get_parameter("flat_cell").as_double();
    prm.viewpoint = {0, 0, 0};  // cloud is in the sensor frame -> sensor origin
    const double range = get_parameter("deviation_range").as_double();

    const planeperc::PipelineOutput out = planeperc::runPipeline(raw, prm);

    std::string summary;
    char buf[320];
    for (size_t k = 0; k < out.planes.size(); ++k) {
      const auto& r = out.planes[k];
      std::snprintf(buf, sizeof(buf),
                    "\n  plane %zu: n=(%.3f %.3f %.3f) d=%.3f m | inliers=%zu (%.0f%%) | flat rms=%.1f mm, peak=%+.1f mm, valley=%+.1f mm",
                    k, r.plane.n.x, r.plane.n.y, r.plane.n.z, r.plane.d, r.inliers, r.inlierRatio * 100,
                    r.flat.rms * 1000, r.flat.maxPos * 1000, r.flat.maxNeg * 1000);
      summary += buf;
    }
    for (const auto& o : planeperc::parallelOffsets(out)) {
      std::snprintf(buf, sizeof(buf),
                    "\n  parallel planes %zu & %zu (%.1f deg apart): plane %zu is %+.1f mm vs plane %zu (+ = closer to sensor)",
                    o.i, o.j, o.angleDeg, o.j, o.offsetM * 1000, o.i);
      summary += buf;
    }
    RCLCPP_INFO_THROTTLE(get_logger(), *get_clock(), 1000, "%zu planes found%s", out.planes.size(), summary.c_str());
    RCLCPP_INFO_THROTTLE(get_logger(), *get_clock(), 1000,
                         "in=%zu voxel=%zu clean=%zu | voxel %.1f SOR %.1f RANSAC %.1f flat %.1f total %.1f ms", out.nIn,
                         out.nVoxel, out.nClean, out.t.voxelMs, out.t.sorMs, out.t.ransacMs, out.t.flatMs, out.t.totalMs);

    // Publish XYZRGB cloud of points that belong to a plane band, coloured by deviation.
    sensor_msgs::msg::PointCloud2 pc;
    pc.header = msg->header;
    pc.height = 1;
    pc.is_bigendian = false;
    pc.is_dense = true;
    pc.point_step = 16;
    auto addField = [&](const char* name, uint32_t off, uint8_t type) {
      sensor_msgs::msg::PointField f;
      f.name = name;
      f.offset = off;
      f.datatype = type;
      f.count = 1;
      pc.fields.push_back(f);
    };
    addField("x", 0, sensor_msgs::msg::PointField::FLOAT32);
    addField("y", 4, sensor_msgs::msg::PointField::FLOAT32);
    addField("z", 8, sensor_msgs::msg::PointField::FLOAT32);
    addField("rgb", 12, sensor_msgs::msg::PointField::FLOAT32);

    size_t cnt = 0;
    for (size_t i = 0; i < out.cleaned.size(); ++i)
      if (!std::isnan(out.deviation[i])) ++cnt;
    pc.width = static_cast<uint32_t>(cnt);
    pc.row_step = pc.point_step * pc.width;
    pc.data.resize(static_cast<size_t>(pc.row_step));
    size_t w = 0;
    for (size_t i = 0; i < out.cleaned.size(); ++i) {
      if (std::isnan(out.deviation[i])) continue;
      uint8_t r, g, b;
      planeperc::deviationColor(out.deviation[i], range, r, g, b);
      const float xyz[3] = {static_cast<float>(out.cleaned[i].x), static_cast<float>(out.cleaned[i].y),
                            static_cast<float>(out.cleaned[i].z)};
      const uint32_t rgb = (static_cast<uint32_t>(r) << 16) | (static_cast<uint32_t>(g) << 8) | b;
      uint8_t* dst = pc.data.data() + w * pc.point_step;
      std::memcpy(dst, xyz, 12);
      std::memcpy(dst + 12, &rgb, 4);
      ++w;
    }
    pub_->publish(pc);
  }

  rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr sub_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_;
};

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<PlaneNode>());
  rclcpp::shutdown();
  return 0;
}
