import matplotlib.pyplot as plt
import numpy as np
import os

def RMSE():
    def read_kf(file_path):
        with open(file_path, 'r') as file:
            data = []
            for line in file.readlines():
                values = list(map(float, line.strip().split('\t')))
                # values[:] = [v * 0.001 for v in values[:]]
                # values[0] += 4.5 
                # values[1] += 4.0  
                # values[2] -= 0.1  

                data.append(values) 
            return data 
                
    def read_data_gt(file_path):
        with open(file_path, 'r') as file:
            data = []
            for line in file.readlines():
                values = list(map(float, line.strip().split('\t')))
                values[:] = [v * 0.001 for v in values[:]]
                values[0] += 4.5 
                values[1] += 4.0  
                # values[2] -= 0.4  

                data.append(values)
            return data      
        
    def interpolate_data(data, target_length):
        x = np.arange(len(data))
        x_new = np.linspace(0, len(data) - 1, target_length)
        return np.array([np.interp(x_new, x, [d[i] for d in data]) for i in range(len(data[0]))]).T

    def calculate_rmse(true_data, pred_data):
        return np.sqrt(np.mean((true_data - pred_data) ** 2, axis=0))

    files = ['../config/hw3_gt.txt', '../config/mppi_pose2.txt']


    uwb_data = read_data_gt(files[0])
    ekf_data = read_kf(files[1])

    max_length = max(len(uwb_data), len(ekf_data))

    uwb_data_interp = interpolate_data(uwb_data, max_length)
    ekf_data_interp = interpolate_data(ekf_data, max_length)


    rmse_ekf = calculate_rmse(uwb_data_interp, ekf_data_interp)


    print("EKF RMSE: ", rmse_ekf[:3])

    fig, axs = plt.subplots(3, 1, figsize=(10, 12))

    axs[0].plot(uwb_data_interp[:, 0], color='black', label='GT', linestyle='-')
    axs[0].plot(ekf_data_interp[:, 0], color='red', label='MPPI', linestyle='-')

    axs[1].plot(uwb_data_interp[:, 1], color='black', label='GT', linestyle='-')
    axs[1].plot(ekf_data_interp[:, 1], color='red', label='MPPI', linestyle='-')

    axs[2].plot(uwb_data_interp[:, 2], color='black', label='GT', linestyle='-')
    axs[2].plot(ekf_data_interp[:, 2], color='red', label='MPPI', linestyle='-')    

    
    axs[0].set_ylabel('X[m]')
    axs[1].set_ylabel('Y[m]')
    axs[2].set_ylabel('Z[m]')
    axs[2].set_xlabel('MPPI Estimation')

    for ax in axs:
        ax.legend(loc='upper right') 
        ax.grid()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    RMSE()
