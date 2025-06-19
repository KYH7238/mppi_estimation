import numpy as np
import matplotlib.pyplot as plt

pose_file = '../data/mppi_pose.txt'
gt_file   = '../data/data1_gt.txt'

def read_pose(path):
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            ts = float(parts[0])
            x, y, z = map(float, parts[1:4])
            data.append([ts, x, y, z])
    return data

def read_gt(path):
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            x = float(parts[0]) * 0.001
            y = float(parts[1]) * 0.001
            z = float(parts[2]) * 0.001
            data.append([x, y, z])
    return data

def interpolate_data(data, length):
    x = np.arange(len(data))
    x_new = np.linspace(0, len(data)-1, length)
    return np.array([
        np.interp(x_new, x, [row[i] for row in data])
        for i in range(len(data[0]))
    ]).T

def calculate_rmse(a, b):
    return np.sqrt(np.mean((a - b)**2, axis=0))

pose_data = read_pose(pose_file)
gt_data   = read_gt(gt_file)

L = max(len(pose_data), len(gt_data))
pose_interp = interpolate_data(pose_data, L)[:,1:4]
gt_interp   = interpolate_data(gt_data,   L)

rmse = calculate_rmse(gt_interp, pose_interp)
print(f'RMSE → X: {rmse[0]:.4f}, Y: {rmse[1]:.4f}, Z: {rmse[2]:.4f}')

fig, axs = plt.subplots(3,1,figsize=(10,12))
axs[0].plot(gt_interp[:,0],   color='black', label='GT')
axs[0].plot(pose_interp[:,0], color='red',   label='Pose')
axs[1].plot(gt_interp[:,1],   color='black', label='GT')
axs[1].plot(pose_interp[:,1], color='red',   label='Pose')
axs[2].plot(gt_interp[:,2],   color='black', label='GT')
axs[2].plot(pose_interp[:,2], color='red',   label='Pose')

axs[0].set_ylabel('X [m]')
axs[1].set_ylabel('Y [m]')
axs[2].set_ylabel('Z [m]')
axs[2].set_xlabel('Sample Index')

for ax in axs:
    ax.legend(loc='upper right')
    ax.grid()

plt.tight_layout()
plt.show()
