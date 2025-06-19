#!/usr/bin/env python3
import rosbag

bagfile = '0612_hw3.bag'
uwb_out = 'uwb_data.txt'
imu_out = 'imu_data.txt'

bag = rosbag.Bag(bagfile)
uwb_msgs = [msg for _, msg, _ in bag.read_messages(topics=['/nlink_linktrack_tagframe0'])]
imu_msgs = [msg for _, msg, _ in bag.read_messages(topics=['/imu/data'])]
bag.close()

with open(uwb_out, 'w') as f:
    for msg in uwb_msgs:
        t = msg.system_time / 1000.0
        ranges = msg.dis_arr
        line = [str(t)] + [str(v) for v in ranges]
        f.write(' '.join(line) + '\n')

with open(imu_out, 'w') as f:
    for msg in imu_msgs:
        t = msg.header.stamp.to_sec()
        gyro = msg.angular_velocity
        acc = msg.linear_acceleration
        line = [str(t), str(gyro.x), str(gyro.y), str(gyro.z),
                str(acc.x), str(acc.y), str(acc.z)]
        f.write(' '.join(line) + '\n')

