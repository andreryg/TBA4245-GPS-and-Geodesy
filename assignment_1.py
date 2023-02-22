"""
Task 1:
Points: Moholt, TP342, ST46
"""
import numpy as np

#Baseline ST46 - TP342:
dx1 = 868.806
dy1 = -2585.885
dz1 = -313.646

#Baseline ST46 - Moholt:
dx2 = 2265.69
dy2 = 248.405
dz2 = -1116.088

#MOHOLT:
latitud = 63.40895119
longitud = 10.43111038

#ST46
latitude = 63.43166874
longitude = 10.43443358

#TP342
latitu = 63.42730197
longitu = 10.38034967



def transformation_matrix(latitude, longitude, dx, dy, dz):
    """
    Calculates delta x, delta y and delta z in a local geodetic system based on a cartesian baseline and the latitude and longitude of the starting point.
    :param latitude: float
    :param longitude: float
    :param dx: float
    :param dy: float
    :param dz: float
    :return: numpy array
    """
    point = np.array([dx,dy,dz]).T
    Trans = np.array([[-np.sin(latitude)*np.cos(longitude), -np.sin(latitude)*np.sin(longitude), np.cos(latitude)],
                      [-np.sin(longitude), np.cos(longitude), 0], 
                        [np.cos(latitude)*np.cos(longitude), np.cos(latitude)*np.sin(longitude), np.sin(latitude)]])
    return np.dot(Trans,point)

def slope_distance(dx, dy, dz):
    return np.sqrt(dx**2+dy**2+dz**2)

def horizontal_distane(dx, dy):
    return np.sqrt(dx**2+dy**2)

def azimuth(dx, dy):
    return np.arctan(dy/dx)
    

local_dx1, local_dy1, local_dz1 = transformation_matrix(latitude, longitude, dx1, dy1, dz1)
local_dx2, local_dy2, local_dz2 = transformation_matrix(latitude, longitude, dx2, dy2, dz2)
print(local_dz1, local_dz2)
#print(transformation_matrix(latitude, longitude, dx2, dy2, dz2))
print(horizontal_distane(local_dx1, local_dy1))
print(azimuth(local_dx1, local_dy1))

