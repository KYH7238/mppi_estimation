import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

def quat_to_euler_deg(qx, qy, qz, qw):
    r = R.from_quat([qx, qy, qz, qw])
    return r.as_euler('xyz', degrees=True)

def rotmat_to_euler_deg(mat):
    mat = np.array(mat).reshape(3,3)
    r = R.from_matrix(mat)
    return r.as_euler('xyz', degrees=True)

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

mppi_lines = open('../config/mppi_pose1.txt','r').readlines()
gt_lines = open('../config/hw3_gt_rotation.txt','r').readlines()

mppi_euler = []
gt_euler = []

for line in mppi_lines:
    vals = list(map(float, line.strip().split()))
    _, _, _, qx, qy, qz, qw = vals
    e_deg = quat_to_euler_deg(qx, qy, qz, qw)
    mppi_euler.append(e_deg)

for line in gt_lines:
    vals = list(map(float, line.strip().split()))
    e_deg = rotmat_to_euler_deg(vals)
    gt_euler.append(e_deg)

mppi_euler = np.array(mppi_euler)
gt_euler = np.array(gt_euler)
max_len = max(len(mppi_euler), len(gt_euler))

mppi_intp = interpolate_data(mppi_euler, max_len)
gt_intp = interpolate_data(gt_euler, max_len)
r = rmse(mppi_intp, gt_intp)

print(f"RMSE(deg): Roll={r[0]:.3f}, Pitch={r[1]:.3f}, Yaw={r[2]:.3f}")

titles = ['Roll(deg)', 'Pitch(deg)', 'Yaw(deg)']
fig, axs = plt.subplots(1, 3, figsize=(12,4))
for i in range(3):
    axs[i].plot(mppi_intp[:, i], label='mppi', color='r')
    axs[i].plot(gt_intp[:, i], label='GT', color='k')
    axs[i].set_title(titles[i])
    axs[i].legend()
    axs[i].grid()
plt.tight_layout()
plt.show()