import numpy as np
import csv

def readfile():
    #empty arrays for storing values
    el_values = []
    m_values = []
    alpha_values = []
    #take file name from user input
    file_name = input("file name to plot: ")

    #open file
    file = open(file_name, "r")
    reader = csv.reader(file)

    #skip header
    next(reader)
    
    #read in format "el,m"
    for row in reader:
        el, m, alpha = map(float, row)
        el_values.append(el)
        m_values.append(m)
        alpha_values.append(alpha)
        

    #sort data in ascending order based on el value, then m value
    data = list(zip(el_values, m_values, alpha_values))
    ordered_data = sorted(data)
    el_values_ordered, m_values_ordered, alpha_values_ordered = zip(*ordered_data)
    #print for test
    print("el:", el_values_ordered)
    print("m:", m_values_ordered)
    print("alpha:", alpha_values_ordered)

    #close file
    file.close()

    return el_values_ordered, m_values_ordered, alpha_values_ordered
