import numpy as np
import matplotlib.pyplot as plt

# Read file from user input
def readfile_Fp(file_name):
    try:
        with open(file_name, "r") as file:
            # Skip the header line
            next(file)
            # Initialize lists to store data
            theta_list = []
            phip_list = []
            Fp_x_list = []
            Fp_y_list = []
            Fp_z_list = []
            # Read each line in the file
            for line in file:
                # Split the line by whitespace
                values = line.split()
                # Extract theta, phip, Fp_x, Fp_y, and Fp_z values
                theta = int(values[0])
                phip = int(values[1])
                Fp_x = float(values[2])
                Fp_y = float(values[3])
                Fp_z = float(values[4])
                # Append values to respective lists
                theta_list.append(theta)
                phip_list.append(phip)
                Fp_x_list.append(np.abs(Fp_x))
                Fp_y_list.append(np.abs(Fp_y))
                Fp_z_list.append(np.abs(Fp_z))
            
            # Convert lists to NumPy arrays
            theta_array = np.array(theta_list)
            phip_array = np.array(phip_list)
            Fp_x_array = np.array(Fp_x_list)
            Fp_y_array = np.array(Fp_y_list)
            Fp_z_array = np.array(Fp_z_list)
            
            # Create a matrix with theta, phip, Fp_x, Fp_y, and Fp_z as columns
            Fp_matrix = np.column_stack((phip_array, theta_array, Fp_x_array, Fp_y_array, Fp_z_array))
            
            return Fp_matrix
    except FileNotFoundError:
        print("File not found")
        return np.array([])


file_name = input("Enter file name: ")
Fp_matrix = readfile_Fp(file_name)


unique_phip_values = np.unique(Fp_matrix[:, 0])


for phip_value in unique_phip_values:

    rows_with_phip_value = Fp_matrix[Fp_matrix[:, 0] == phip_value]

    theta_values = rows_with_phip_value[:, 1]
    Fp_x_values = rows_with_phip_value[:, 2]
    Fp_y_values = rows_with_phip_value[:, 3]
    Fp_z_values = rows_with_phip_value[:, 4]

    # Plot
    plt.figure()
    plt.plot(theta_values, Fp_x_values, marker='o', color='blue', label='Fp_x')
    plt.plot(theta_values, Fp_y_values, marker='o', color='green', label='Fp_y')
    plt.plot(theta_values, Fp_z_values, marker='o', color='red', label='Fp_z')
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'${Fp}$')
    plt.title('Fp graph for psi = {}'.format(phip_value))
    plt.legend()
    plt.grid(True)
    plt.show()
