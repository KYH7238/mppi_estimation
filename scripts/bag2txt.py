#!/usr/bin/env python3
import sys
import rosbag
from sensor_msgs.msg import Imu
from nlink_parser.msg import LinktrackTagframe0

def dump_bag(bag_path, imu_out, uwb_out):
    bag = rosbag.Bag(bag_path, 'r')

    with open(imu_out, 'w') as f_imu:
        for topic, msg, t in bag.read_messages(topics=['/vectornav_driver_node/imu/data']):
            ts = msg.header.stamp.secs + msg.header.stamp.nsecs * 1e-9
            ax = msg.linear_acceleration.x
            ay = msg.linear_acceleration.y
            az = msg.linear_acceleration.z
            gx = msg.angular_velocity.x
            gy = msg.angular_velocity.y
            gz = msg.angular_velocity.z
            f_imu.write(f"{ts:.6f} {ax:.6f} {ay:.6f} {az:.6f} "
                        f"{gx:.6f} {gy:.6f} {gz:.6f}\n")

    with open(uwb_out, 'w') as f_uwb:
        for topic, msg, t in bag.read_messages(topics=['/nlink_linktrack_tagframe0']):
            ts = msg.system_time / 1000.0
            tag_id = msg.id
            ranges = msg.dis_arr
            ranges_str = ' '.join(f"{r:.6f}" for r in ranges)
            f_uwb.write(f"{tag_id} {ts:.6f} {ranges_str}\n")

    bag.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.bag> <imu.txt> <uwb.txt>")
        sys.exit(1)
    bag_path = sys.argv[1]
    imu_txt   = sys.argv[2]
    uwb_txt   = sys.argv[3]
    dump_bag(bag_path, imu_txt, uwb_txt)
    print(f"Done: IMU → {imu_txt}, UWB → {uwb_txt}")
