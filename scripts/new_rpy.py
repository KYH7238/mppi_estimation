import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

def read_gt_rpy(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            vals = list(map(float, line.strip().split()))
            rpy = vals[3:6]
            data.append(rpy)
    return np.array(data)

def quat_to_euler_deg(qx, qy, qz, qw):
    r = R.from_quat([qx, qy, qz, qw])
    return r.as_euler('xyz', degrees=True)

def read_mppi_rpy(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            vals = list(map(float, line.strip().split()))
            qx, qy, qz, qw = vals[3:7]
            rpy_deg = quat_to_euler_deg(qx, qy, qz, qw)
            data.append(rpy_deg)
    return np.array(data)

def interpolate_data(data, target_length):
    data = np.array(data)
    if data.shape[0] < 2:
        return data
    x_old = np.arange(data.shape[0])
    x_new = np.linspace(0, data.shape[0]-1, target_length)
    result = []
    for col in range(data.shape[1]):
        y_new = np.interp(x_new, x_old, data[:,col])
        result.append(y_new)
    return np.array(result).T

def rmse(a, b):
    return np.sqrt(np.mean((a - b)**2, axis=0))

gt_file = '../config/20250527/0502_hw1_gt.txt'
mppi_file = '../config/20250527/uwb_node.txt'
# gt_file = '../config/hw3_gt.txt'
# mppi_file = '../config/mppi_pose1.txt'

gt_rpy = read_gt_rpy(gt_file)
mppi_rpy = read_mppi_rpy(mppi_file)

max_len = max(len(gt_rpy), len(mppi_rpy))
gt_rpy_interp = interpolate_data(gt_rpy, max_len)
mppi_rpy_interp = interpolate_data(mppi_rpy, max_len)

r = rmse(mppi_rpy_interp, gt_rpy_interp)
print(f"RMSE(deg): Roll={r[0]:.3f}, Pitch={r[1]:.3f}, Yaw={r[2]:.3f}")

titles = ['Roll(deg)', 'Pitch(deg)', 'Yaw(deg)']
fig, axs = plt.subplots(1, 3, figsize=(12,4))
for i in range(3):
    axs[i].plot(mppi_rpy_interp[:, i], label='mppi', color='r')
    axs[i].plot(gt_rpy_interp[:, i], label='GT', color='k')
    axs[i].set_title(titles[i])
    axs[i].legend()
    axs[i].grid()
plt.tight_layout()
plt.show() 