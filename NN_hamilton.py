#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 14:08:54 2022

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

def NearestNeighbour_hamilton(C):
    
    """ Compute Hamiltonian path from input coordinates according to Nearest 
        Neighbour algorithm
    
    Input:
    C -- nested list with x and y coordinates (list)
    
    Output:
    Q -- ordered list of visited nodes (list)
    W -- path length in km (float)
    """
    
    # Initialize used parameter S (list of nodes' status, W, Q)
    S = ["N"] * (len(C)+1)     
    W = 0                                
    Q = []       

    # Select starting node               
    u = choice(C)
    u = C.index(u)
    
    # Save initial node u(s)
    int_u = u
    
    # Loop through input list C until all nodes are PROCCESED ("P")    
    while S[u]=="N":  
        
        # Change initial/current node's status and append it to Q
        S[u] = "P"
        Q.append(u) 
        
        # Extract node's coordinates
        x_u, y_u = C[u]
        
        # Initialize change of W 
        dW = float("Inf")
    
        # Loop through input list C to find UNPROCCESED ("N") nodes  
        for i, el in enumerate(C):
            if S[i] == 'P':
                continue
            else:
                
                # Extract coordinates of new node v(i)  
                x_v, y_v = el
                
                # Compute distance between initial/current node u(s) and new node v(i)  
                d = sqrt((x_u-x_v)**2 + (y_u-y_v)**2)
                
                # Find smallest change W (dW) and set node v(i) as starting node u(s)
                if d < dW:
                    dW = d
                    u = i
                    
        # Add dW to W                
        if dW != float("Inf"):
            W += dW
        
    # Compute and add distance between last and initial node to W      
    last_d = sqrt((C[int_u][0]-C[u][0])**2 + (C[int_u][1]-C[u][1])**2)
    W += last_d
    
    # Change W unit from m to km
    W = W/1000
    
    # Return Output 
    return Q, W


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

def iteration_NearestNeighbour(C, max_iter):
    
    """ Iter Nearest Neighbour algorithm i-times to find best result
        and return minimal path length
    
    Input:
    C -- nested list with x and y coordinates (list)
    max_iter -- nmaximal number of iterations
    
    Return:
    D_results -- list of computed paths' lengths (list)
    Q_optimal -- index of initial node for path with min lenght W_min
    W_min -- minimal path length in m (float)
    """
    # Initialize Output variables
    D_results = []
    Q_initial = []
    
    # Repeat NN algorithm i-times and add results to lists
    for i in range(max_iter):
        
        Q, W = NearestNeighbour_hamilton(C)
        D_results.append(W)
        Q_initial.append(Q)
    
    # Find minimal path length 
    W_min = min(D_results)
    Q_optimal = Q_initial[D_results.index(W_min)][0]
    
    # Return Outputs
    return D_results, Q_optimal, W_min     


# Call all functions with desired Inputs 
# files -- 'restaurants_coord.txt' or 'bus_coordinates'
coordinates = openTXT_exportCoordinates('restaurants_coord.txt')
Q, W = NearestNeighbour_hamilton(coordinates)
print("Path length W is around: {:.3f} km".format(W))
P = path(Q, coordinates)
visualize_hamiltionianP(coordinates, P) 

# # Save starting time before itering 
# start = time.perf_counter()
     
# D_results, Q_optimal, W_min = iteration_NearestNeighbour(coordinates, 200) 

# # Save termination time
# stop = time.perf_counter()

# # Measure complete duration of the code execution
# duration = stop - start

# print("Minimal computed path length W for i random iteration is around: {:.3f} km".format(W_min))
# print("Initial node for best found Hamiltonian path is: {} ".format(Q_optimal))
# print("Duration of the code execution is: {} s".format(duration))





