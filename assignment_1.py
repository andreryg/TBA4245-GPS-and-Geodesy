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

test1 = 63.40919533
test2 = 10.58311814

#test_dx = 

#ST46-TP342 correlation matrix
std_st46_tp342 = [0.0002, 0.0001, 0.0004]
K_st46_tp342 = np.array([[1.0000, 0.2226, 0.2668],
                         [0.2226, 1.0000, 0.2386], 
                         [0.2668, 0.2386, 1.0000]])


def transformation_matrix(latitude, longitude):
    """
    Calculates the transformation matrix required to transform a point defined by latitude and longitude to a local geodetic system.
    :param latitude: float
    :param longitude: float
    :return: 2d numpy array
    """
    latitude = np.deg2rad(latitude)
    longitude = np.deg2rad(longitude)
    Trans = np.array([[-np.sin(latitude)*np.cos(longitude), -np.sin(latitude)*np.sin(longitude), np.cos(latitude)],
                      [-np.sin(longitude), np.cos(longitude), 0], 
                        [np.cos(latitude)*np.cos(longitude), np.cos(latitude)*np.sin(longitude), np.sin(latitude)]])
    print(Trans)
    return Trans

def transform(T, dx, dy, dz):
    """
    Transform a baseline to a local geodetic system.
    :param T: transformation matrix
    :param dx: float
    :param dy: float
    :param dz: float
    :return: numpy array
    """
    point = np.array([[dx],[dy],[dz]])
    return np.dot(T, point)

def slope_distance(dx, dy, dz):
    return np.sqrt(dx**2+dy**2+dz**2)

def horizontal_distane(dx, dy):
    return np.sqrt(dx**2+dy**2)

def azimuth(dx, dy):
    return np.rad2deg(np.arctan(dy/dx))

def covariance_matrix(std, K):
    std_matrix = np.array([[std[0], 0, 0],
                             [0, std[1], 0], 
                             [0, 0, std[2]]])
    return np.dot(np.dot(std_matrix, K), std_matrix)

def distance_UTM(Shh, R, Hs, z):
    return np.arctan(Shh / (R + Hs + z))
    

if __name__ == "__main__":
    local_dx1, local_dy1, local_dz1 = transform(transformation_matrix(latitude, longitude), dx1, dy1, dz1)
    local_dx2, local_dy2, local_dz2 = transform(transformation_matrix(latitude, longitude), dx2, dy2, dz2)
    print(local_dx1, local_dy1)
    print(local_dz1, local_dz2)

    #print(transformation_matrix(latitude, longitude, dx2, dy2, dz2))
    #print(distance_UTM((horizontal_distane(local_dx1, local_dy1)), , )
    print(azimuth(local_dx1, local_dy1))
    print(azimuth(local_dx2, local_dy2))


    transform(transformation_matrix(test1, test2), dx1, dy1, dz1)

    print(covariance_matrix(std_st46_tp342, K_st46_tp342))
