import numpy as np
import matplotlib.pyplot as plt

def readfile_Tp(file_name):
    try:
        with open(file_name, "r") as file:
            # Initialize empty dictionary to store Fp values for each theta
            Tp_data = {}
            # Skip the header line
            next(file)
            theta_list = []
            phip_list = []
            Tp_x_list = []
            Tp_y_list = []
            Tp_z_list = []
            # Read each line in the file
            for line in file:
                # Split the line by whitespace
                values = line.split()
                # Extract theta, phip, Fp_x, Fp_y, and Fp_z values
                theta = int(values[0])
                phip = int(values[1])
                Tp_x = float(values[5])
                Tp_y = float(values[6])
                Tp_z = float(values[7])
                  # Append values to respective lists
                theta_list.append(theta)
                phip_list.append(phip)
                Tp_x_list.append(np.abs(Tp_x))
                Tp_y_list.append(np.abs(Tp_y))
                Tp_z_list.append(np.abs(Tp_z))  

            theta_array = np.array(theta_list)
            phip_array = np.array(phip_list)
            Tp_x_array = np.array(Tp_x_list)
            Tp_y_array = np.array(Tp_y_list)
            Tp_z_array = np.array(Tp_z_list)
            
            # Create a matrix with theta, phip, Fp_x, Fp_y, and Fp_z as columns
            Tp_matrix = np.column_stack((phip_array, theta_array, Tp_x_array, Tp_y_array,Tp_z_array))
            

            return Tp_matrix
    except FileNotFoundError:
        print("File not found")
        return np.array([])


file_name = input("Enter file name: ")
Tp_matrix = readfile_Tp(file_name)

unique_phip_values = np.unique(Tp_matrix[:, 0])


for phip_value in unique_phip_values:

    rows_with_phip_value = Tp_matrix[Tp_matrix[:, 0] == phip_value]

    theta_values = rows_with_phip_value[:, 1]
    Tp_x_values = rows_with_phip_value[:, 2]
    Tp_y_values = rows_with_phip_value[:, 3]
    Tp_z_values = rows_with_phip_value[:, 4]

    # Plot
    plt.figure()
    plt.plot(theta_values, Tp_x_values, marker='o', color='blue', label='Tp_x')
    plt.plot(theta_values, Tp_y_values, marker='o', color='green', label='Tp_y')
    plt.plot(theta_values, Tp_z_values, marker='o', color='red', label='Tp_z')
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'${Fp}$')
    plt.title('Fp graph for psi = {}'.format(phip_value))
    plt.legend()
    plt.grid(True)
    plt.show()