import rospy
import rosbag

position_x = []
position_y = []
position_z = []
orientation_x = []
orientation_y = []
orientation_z = []
orientation_w = []

def callback_path(msg):
    global position_x, position_y, position_z, orientation_x, orientation_y, orientation_z, orientation_w
    position_x.append(msg.pose.position.x)
    position_y.append(msg.pose.position.y)
    position_z.append(msg.pose.position.z)
    orientation_x.append(msg.pose.orientation.x)
    orientation_y.append(msg.pose.orientation.y)
    orientation_z.append(msg.pose.orientation.z)
    orientation_w.append(msg.pose.orientation.w)

def make_txt_file_for_plot(file_name):
    global position_x, position_y, position_z, orientation_x, orientation_y, orientation_z, orientation_w
    with open(file_name, "w") as file:
        for a, b, c, d, e, f, g in zip(position_x, position_y, position_z, orientation_x, orientation_y, orientation_z, orientation_w):
            file.write(f"{a}\t{b}\t{c}\t{d}\t{e}\t{f}\t{g}\n")

if __name__ == "__main__":
    file_name_1 = "../config/mppi_pose1.txt"

    bag_path = "/home/kim/drone_ws/src/mppi_estimation/bag/2025-02-26-15-45-06.bag"
    topics = [("/mppi_pose", file_name_1)]
    for topic_name, file_name in topics:

        position_x.clear()
        position_y.clear()
        position_z.clear()
        orientation_x.clear()
        orientation_y.clear()
        orientation_z.clear()
        orientation_w.clear()

        with rosbag.Bag(bag_path, "r") as bag:
            for topic, msg, t in bag.read_messages(topics=[topic_name]):
                callback_path(msg)

        make_txt_file_for_plot(file_name)
