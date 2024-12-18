import numpy as np
import matplotlib.pyplot as plt

# Read file from user input
def readfile_Fp(file_name):
    try:
        with open(file_name, "r") as file:
            # Skip the header line
            next(file)
            # Initialize lists to store data
            Ompy_list = []
            Up_list = []
            Fp_z_list = []
            Tp_y_list = []
            # Read each line in the file
            for line in file:
                # Split the line by whitespace
                values = line.split()
                # Extract Ompy, Up, Fp_z, and Tp_y values
                Ompy = float(values[0])
                Up = float(values[1])
                Fp_z = float(values[4])
                Tp_y = float(values[6])
                # Append values to respective lists
                Ompy_list.append(Ompy)
                Up_list.append(Up)
                Fp_z_list.append(abs(Fp_z))
                Tp_y_list.append(abs(Tp_y))
            
            # Convert lists to NumPy arrays
            Ompy_array = np.array(Ompy_list)
            Up_array = np.array(Up_list)
            Fp_z_array = np.array(Fp_z_list)
            Tp_y_array = np.array(Tp_y_list)
            
            # Create a matrix with Ompy, Up, Fp_z, and Tp_y as columns
            Fp_matrix = np.column_stack((Ompy_array, Up_array, Fp_z_array, Tp_y_array))
            
            return Fp_matrix
    except FileNotFoundError:
        print("File not found")
        return np.array([])


file_name = input("Enter file name: ")
Fp_matrix = readfile_Fp(file_name)

unique_ompy_values = np.unique(Fp_matrix[:, 0])

for ompy_value in unique_ompy_values:
    # Filter data for each unique Ompy value
    filtered_data = Fp_matrix[Fp_matrix[:, 0] == ompy_value]

    Up_values = filtered_data[:, 1]
    Fp_z_values = filtered_data[:, 2]
    Tp_y_values = filtered_data[:, 3]

    # Plot
    plt.figure()
    plt.plot(Up_values, Fp_z_values, marker='o', color='blue', label='Fp_z')
    plt.plot(Up_values, Tp_y_values, marker='o', color='green', label='Tp_y')
    plt.xlabel('Up')
    plt.ylabel('Value')
    plt.title('Fp_z and Tp_y as a function of Up (Ompy = {})'.format(ompy_value))
    plt.legend()
    plt.grid(True)
    plt.show()
