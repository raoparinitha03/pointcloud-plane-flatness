#!/usr/bin/env python3
"""Publishes a synthetic noisy room (floor + 2 walls + a 30 mm bump on one wall + outliers)
as sensor_msgs/PointCloud2 on /points, with the 'sensor' at the room centre.
Use it to test plane_node without a depth camera in Gazebo.
Run:  python3 fake_room_publisher.py
"""
import numpy as np
import rclpy
from rclpy.node import Node
from sensor_msgs.msg import PointCloud2, PointField

X, Y, Z = 5.0, 4.0, 2.6
DENSITY = 3000        # points per m^2
NOISE = 0.004         # 4 mm depth noise
BUMP_AMP, BUMP_SIGMA, BUMP_Y, BUMP_Z = 0.03, 0.15, 2.0, 1.3


def make_room(rng):
    n = lambda area: int(area * DENSITY)
    nf, na, nb = n(X * Y), n(Y * Z), n(X * Z)
    floor = np.c_[rng.uniform(0, X, nf), rng.uniform(0, Y, nf), rng.normal(0, NOISE, nf)]
    ya, za = rng.uniform(0, Y, na), rng.uniform(0, Z, na)
    r2 = (ya - BUMP_Y) ** 2 + (za - BUMP_Z) ** 2
    wall_a = np.c_[rng.normal(0, NOISE, na) + BUMP_AMP * np.exp(-r2 / (2 * BUMP_SIGMA ** 2)), ya, za]
    wall_b = np.c_[rng.uniform(0, X, nb), rng.normal(0, NOISE, nb), rng.uniform(0, Z, nb)]
    pts = np.vstack([floor, wall_a, wall_b])
    no = int(0.03 * len(pts))
    junk = np.c_[rng.uniform(-0.5, X + 0.5, no), rng.uniform(-0.5, Y + 0.5, no), rng.uniform(-0.2, Z + 0.2, no)]
    pts = np.vstack([pts, junk]) - np.array([X / 2, Y / 2, Z / 2])  # sensor at room centre
    return pts.astype(np.float32)


class Fake(Node):
    def __init__(self):
        super().__init__('fake_room_publisher')
        self.pub = self.create_publisher(PointCloud2, 'points', 1)
        self.rng = np.random.default_rng(0)
        self.create_timer(0.5, self.tick)

    def tick(self):
        pts = make_room(self.rng)
        msg = PointCloud2()
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.header.frame_id = 'map'
        msg.height, msg.width = 1, len(pts)
        msg.fields = [PointField(name=c, offset=4 * i, datatype=PointField.FLOAT32, count=1) for i, c in enumerate('xyz')]
        msg.is_bigendian, msg.point_step, msg.is_dense = False, 12, True
        msg.row_step = 12 * len(pts)
        msg.data = pts.tobytes()
        self.pub.publish(msg)


def main():
    rclpy.init()
    rclpy.spin(Fake())


if __name__ == '__main__':
    main()
