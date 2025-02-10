#!/usr/bin/env python3
import rospy
from nav_msgs.msg import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

trajectory_data = []

def path_callback(msg):
    global trajectory_data
    trajectory_data = []
    for pose_stamped in msg.poses:
        x = pose_stamped.pose.position.x
        y = pose_stamped.pose.position.y
        trajectory_data.append((x, y))

def animate(frame):
    plt.cla()  
    if trajectory_data:
        xs, ys = zip(*trajectory_data)
        plt.plot(xs, ys, 'b-', marker='o')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.xlim(0,20)
        plt.ylim(0,20)
        plt.title("Trajectory")
    else:
        plt.text(0.5, 0.5, "No Data", horizontalalignment='center')

def main():
    rospy.init_node('trajectory_plotter', anonymous=True)
    rospy.Subscriber("trajectory", Path, path_callback)

    ani = FuncAnimation(plt.gcf(), animate, interval=100)
    plt.show()

if __name__ == '__main__':
    main()
