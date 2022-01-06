#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 10:23:53 2022

@author: juliasuskova
"""

import math
from math import *
import random 
from random import *
import numpy as np
import matplotlib.pyplot as plt
import csv
import time


def openTXT_exportCoordinates(filename):
    
    """ Open txt file as csv and export list of coordinates
    
    Input:
    filename -- file path (string)
    
    Output:
    coordinates -- list with x and y coordinates (list)
    """
    
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=";")
        coordinates = []
        for row in csv_reader:
            row[0] = float(row[0])
            row[1] = float(row[1])
            coordinates.append(row)
            
    return coordinates

def BestInsertion_hamilton(C):
    
    """ Compute Hamiltonian path from input coordinates according to Best
        Insertion algorithm
    
    Input:
    C -- nested list with x and y coordinates (list)
    
    Output:
    W -- path length in km (float)
    modified_path -- modified path of visited nodes (list)
    """
    
    # Extract x and y from coordinates nested list
    x = [x[0] for x in C]
    y = [y[1] for y in C]
    
    # Choose 3 random initial nodes  
    indices = range(len(coordinates))
    init_path = sample(indices, 3)
    
    
    # Close initial path
    init_path.append(init_path[0])
    
    # Initialize path length W
    W = 0
    for i in range(len(init_path)-1):
        dW = sqrt((x[init_path[i]] - x[init_path[i + 1]])**2 + (y[init_path[i]] - y[init_path[i + 1]])**2)    
        W += dW
    
    # Remove initial path's nodes from list C
    C_copy = C.copy()
    for i in range(3):
        C_copy.remove(C[init_path[i]])
    
    # Initialize node u from C_copy and extract its index
    u = choice(C_copy)
    u = C.index(u)
    
    # Initialize used parameter S (list of nodes' status) 
    S = ["N"] * (len(C)+1)
    
    # Loop through copy list C_copy until all nodes are PROCCESED ("P")  
    while S[u] == "N":
        S[u] = "P"
        
        # Initialize change of path length 
        dW = float("Inf")

        # Loop through initial path and extract disatnces
        for i in range(len(init_path)-1):
            
            # Initialize node s and node v
            s = init_path[i]
            v = init_path[i+1]
            
            # Calculate dW(s,u), dW(u,v) and dW(s,v)
            dW_su = sqrt((x[s] - x[u])**2 + (y[s] - y[u])**2)
            dW_uv = sqrt((x[v] - x[u])**2 + (y[v] - y[u])**2)            
            dW_sv = sqrt((x[v] - x[s])**2 + (y[v] - y[s])**2)

            # Calculate potential updated W
            W_potential = W + dW_su + dW_uv - dW_sv
            if W_potential < dW:
                dW = W_potential
                index = i+1
        
        # Update W
        W = dW
        
        # Add new node u to the path in the correct order
        init_path.insert(index, u)
        
        # Remove node u from C_copy and choose new u until len(C_copy) == 1 
        C_copy.remove(C[u])
        
        # Select next random node u until only one node ist left
        if len(C_copy) > 1:
            u = choice(C_copy)  
            u = C.index(u)
        else:  
            u = C.index(C_copy[0])
            S[u] = "P"
    

    # Initialize change of path length 
    dW = float("Inf")
    
    # Find best positioning for the last node u
    for i in range(len(init_path)-1):     
        
        # Initialize node s and node v
        s = init_path[i]
        v = init_path[i+1]
        
        # Calculate dW(s,u), dW(u,v) and dW(s,v)
        dW_su = sqrt((x[s] - x[u])**2 + (y[s] - y[u])**2)
        dW_uv = sqrt((x[v] - x[u])**2 + (y[v] - y[u])**2)            
        dW_sv = sqrt((x[v] - x[s])**2 + (y[v] - y[s])**2)

        # Calculate potential updated W
        W_potential = W + dW_su + dW_uv - dW_sv
        if W_potential < dW:
            dW = W_potential
            index = i+1
    
    # Update W and change from [m] to [km]
    W = dW
    W = W/1000
        
    # Add new node u to the path in the correct order
    init_path.insert(index, u)
    modified_path = init_path
    
    
    return W, modified_path


def path(Q, C):
        
    """ Extract coordinates of visited nodes in correct order for later 
        visualization
    
    Input:
    Q -- ordered list of visited nodes (list)
    C -- nested list with x and y coordinates (list)
    
    Output:
    P -- ordered list of nodes' coordinates (list)
    """
    
    # Initialize P
    P = []

    # Loop through Q to extract coordinates    
    for i in Q:
        P.append(C[i])
    
    # Add initial node's coordinates to P
    P.append(C[Q[0]])
    
    # Change list P to array P to visualize edges
    P = np.array(P)

    # Return Output
    return P

def visualize_hamiltionianP(C, P):
    
    """ Extract coordinates of visited nodes in correct order for later 
        visualization
    
    Input:
    C -- nested list with x and y coordinates (list)
    P -- ordered list of nodes' coordinates (list)
    
    Return:
    graphical visualization of computed Hamiltonian path
    """
    
    # Store x and y coordinates into separate variables 
    x = [x[0] for x in C]
    y = [y[1] for y in C]
    
    # Display the graph
    plt.scatter(x, y, c='red')
    plt.plot(P[:, 0], P[:, 1], c='black')
    plt.show()
    
def iteration_BestInsertion(C, max_iter):
    
    """ Iter Best Insertion algorithm i-times to find best result
        and return minimal path length
    
    Input:
    C -- nested list with x and y coordinates (list)
    max_iter -- nmaximal number of iterations
    
    Return:
    D_results -- list of computed paths' lengths (list)
    W_min -- minimal path length in m (float)
    P_optimal -- optimal path with minimal found lenght (list)  
    """
    # Initialize Output variables
    D_results = []
    Paths = []
    
    # Repeat BI algorithm i-times and add results to lists
    for i in range(max_iter):
        
        # Apply BestInsertion algorithm
        W, P = BestInsertion_hamilton(C)
        D_results.append(W)
        Paths.append(P)
    
    # Find minimal path length 
    W_min = min(D_results)
    P_optimal = Paths[D_results.index(W_min)]
    
    
    # Return Outputs
    return D_results, W_min, P_optimal  


# Call all functions with desired Inputs 
# files -- 'restaurants_coord.txt' or 'bus_coordinates'
coordinates = openTXT_exportCoordinates('restaurants_coord.txt')
W, P = BestInsertion_hamilton(coordinates)
print("Path length W is around: {:.3f} km".format(W))
P = path(P, coordinates)
visualize_hamiltionianP(coordinates, P) 

# # Save starting time before itering 
# start = time.perf_counter()

  
# D_results, W_min, P_optimal = iteration_BestInsertion(coordinates, 200) 

# # Save termination time
# stop = time.perf_counter()

# # Measure complete duration of the code execution
# duration = stop - start

# print("Minimal computed path length W for i random iterations is around: {:.3f} km".format(W_min))
# print("Duration of the code execution is: {} s".format(duration))



