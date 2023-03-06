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

a = 6378137
b = 6356752.3141

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
    #print(Trans)
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
    """
    Calculates azimuth in gradian, angle from north, based on a baseline.
    """
    if dx > 0 and dy > 0: #Zone 1
        return np.arctan(dy/dx)*200/np.pi 
    elif dx <= 0 and dy > 0: #Zone 2
        return -np.arctan(dy/dx)*200/np.pi + 200 
    elif dx <= 0 and dy <= 0: #Zone 3
        return np.arctan(dy/dx)*200/np.pi + 200 
    else: #Zone 4
        return -np.arctan(dy/dx)*200/np.pi + 400

def covariance_matrix(std, K):
    std_matrix = np.array([[std[0], 0, 0],
                             [0, std[1], 0], 
                             [0, 0, std[2]]])
    return np.dot(np.dot(std_matrix, K), std_matrix)

def distance_UTM(Shh, R, Hs, z, ya, yb):
    S0 = R * np.arctan(Shh / (R + Hs + z))
    Sm = S0 + S0/(6*R**2) * (ya**2 + yb**2 + ya * yb) - 0.0004 * S0
    
    return Sm

def radius(latitude, A, a = 6378137, b = 6356752.3141):
    e_2 = 1 - b**2/a**2
    w = (1-e_2*np.sin(np.deg2rad(latitude))**2)**(1/2)
    rN = a / w
    rM = a * (1-e_2) / (w ** 3)
    R = (rM * rN) / (rN*np.cos(A*np.pi/200)**2 + rM*np.sin(A*np.pi/200)**2)
    return R

def UTM_bearing(latitude, longitude, R, xa, xb, ya, yb, A, a = 6378137, b = 6356752.3141):
    e_2 = 1 - b**2/a**2
    long = np.deg2rad(longitude-9)
    lat = np.deg2rad(latitude)
    eta_2 = e_2 * np.cos(lat)**2 / (1-e_2)
    c = long * np.sin(lat) + long**3/3 * np.sin(lat) * np.cos(lat)**2 * (1 + 3*eta_2 + 2*eta_2**2) + long**5/5 * np.sin(lat) * np.cos(lat)**4 * (2 - np.tan(lat)**2)

    ro = 200/np.pi
    c = c * ro

    delta = ro / (6 * R**2) * (2*ya+ yb) * (xb - xa)
    #A = A * 400/360
    bearing = A - abs(delta) - abs(c)

    return bearing

def main():
    #TASK 1
    print(f"----------------------Task 1a----------------------")
    local_dx1, local_dy1, local_dz1 = transform(transformation_matrix(latitude, longitude), dx1, dy1, dz1) #TP342
    local_dx2, local_dy2, local_dz2 = transform(transformation_matrix(latitude, longitude), dx2, dy2, dz2) #MOHOLT
    print(f"ST46-TP342: {[local_dx1[0], local_dy1[0], local_dz1[0]]}")
    print(f"ST46-MOHOLT: {[local_dx2[0], local_dy2[0], local_dz2[0]]}")
    

    #print(local_dx1, local_dy1)
    #print(local_dz1, local_dz2)


    #print("r = ", radius(latitude, azimuth(local_dx1, local_dy1)))

    #TASK 1b
    print(f"----------------------Task 1b----------------------")
    utm_dist_1 = distance_UTM(horizontal_distane(local_dx1, local_dy1), radius(latitude, azimuth(local_dx1, local_dy1)), slope_distance(local_dx1, local_dy1, local_dz1), local_dz1, 0, local_dy1)
    utm_dist_2 = distance_UTM(horizontal_distane(local_dx2, local_dy2), radius(latitude, azimuth(local_dx2, local_dy2)), slope_distance(local_dx2, local_dy2, local_dz2), local_dz2, 0, local_dy2)
    print(f"ST46-TP342 UTM distance: {utm_dist_1[0]} m, compared to values from table 2: {np.sqrt((7033941.628-7034487.402)**2 + (568890.318 - 571578.304)**2)} m, Difference: {utm_dist_1[0] - np.sqrt((7033941.628-7034487.402)**2 + (568890.318 - 571578.304)**2)} m.")
    print(f"ST46-MOHOLT UTM distance: {utm_dist_2[0]} m, compared to values from table 2: {np.sqrt((7031952.892-7034487.402)**2 + (571469.041 - 571578.304)**2)} m, Difference: {utm_dist_2[0] - np.sqrt((7031952.892-7034487.402)**2 + (571469.041 - 571578.304)**2)} m.")

    #TASK 1c
    print(f"----------------------Task 1c----------------------")
    bearing_1 = UTM_bearing(latitude, longitude, radius(latitude, azimuth(local_dx1, local_dy1)), 0, local_dx1, 0, local_dy1, azimuth(local_dx1, local_dy1))
    bearing_2 = UTM_bearing(latitude, longitude, radius(latitude, azimuth(local_dx2, local_dy2)), 0, local_dx2, 0, local_dy2, azimuth(local_dx2, local_dy2))
    print(f"ST46-TP342 bearing: {bearing_1[0]} grad, compared to values from table 2: {azimuth(7033941.628-7034487.402, 568890.318 - 571578.304)} grad, Difference: {bearing_1[0] - azimuth(7033941.628-7034487.402, 568890.318 - 571578.304)} grad.")
    print(f"ST46-MOHOLT bearing: {bearing_2[0]} grad, compared to values from table 2: {azimuth(7031952.892-7034487.402, 571469.041 - 571578.304)} grad, Difference: {bearing_2[0] - azimuth(7031952.892-7034487.402, 571469.041 - 571578.304)} grad.")

    #TASK 1d
    print(f"----------------------Task 1c----------------------")
    #print(transformation_matrix(latitude, longitude, dx2, dy2, dz2))
    #print(distance_UTM((horizontal_distane(local_dx1, local_dy1)), , )
    #print(azimuth(local_dx1, local_dy1))
    #print(azimuth(local_dx2, local_dy2))


    #transform(transformation_matrix(test1, test2), dx1, dy1, dz1)

    #print(covariance_matrix(std_st46_tp342, K_st46_tp342))


if __name__ == "__main__":
    main()