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

N_moholt = 170.839 - 131.404
N_tp342 = 44.618 - 5.194
N_st46 = 151.851 - 112.554

#test_dx = 

#ST46-TP342 correlation matrix
std_st46_tp342 = [0.0002, 0.0001, 0.0004]
K_st46_tp342 = np.array([[1.0000, 0.2226, 0.2668],
                         [0.2226, 1.0000, 0.2386], 
                         [0.2668, 0.2386, 1.0000]])
std_st46_moholt = [0.0001, 0.0001, 0.0003]
K_std46_moholt = np.array([[1.0000, 0.3026, 0.3727],
                           [0.3026, 1.0000, 0.2288],
                           [0.3727, 0.2288, 1.0000]])
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
    
def zenith(dx, dy, dz):
    """
    Calculates zenith angle in gradian between two points.
    """
    return np.arctan(horizontal_distane(dx,dy)/dz)*200/np.pi

def covariance_matrix(std, K):
    """
    Calculates variance/covariance matrix [mm2] from standard deviations and correlation matrix of baseline. 
    """
    std_matrix = np.array([[std[0]*1000, 0, 0],
                             [0, std[1]*1000, 0], 
                             [0, 0, std[2]*1000]])
    return std_matrix @ K @ std_matrix

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

def NN2000_height_diff(zv, dx, dy, dz, N1, N2, R):
    dN = N2 - N1
    dHE = (slope_distance(dx, dy, dz) * np.sin(zv*np.pi/200))**2 / (2*R)

    dH = dz - dN + dHE
    return dH

def error_prop_law(A, B):
    return A @ B @ A.T

def F_matrix(A, Sh):
    A = A * np.pi/200
    Sh = Sh * 1000
    rho = 200000/np.pi
    F11 = -np.sin(A)*rho / Sh
    F12 = np.cos(A)*rho / Sh
    F21 = np.cos(A)
    F22 = np.sin(A)
    F = np.array([[float(F11), float(F12), 0.],
                  [float(F21), float(F22), 0.],
                  [0., 0., 1.]])
    return F


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
    print(f"ST46-TP342 UTM distance: {utm_dist_1[0]} m, compared to values from table 2: {np.sqrt((7033941.628-7034487.402)**2 + (568890.318 - 571578.304)**2)} m. Difference: {utm_dist_1[0] - np.sqrt((7033941.628-7034487.402)**2 + (568890.318 - 571578.304)**2)} m.")
    print(f"ST46-MOHOLT UTM distance: {utm_dist_2[0]} m, compared to values from table 2: {np.sqrt((7031952.892-7034487.402)**2 + (571469.041 - 571578.304)**2)} m. Difference: {utm_dist_2[0] - np.sqrt((7031952.892-7034487.402)**2 + (571469.041 - 571578.304)**2)} m.")

    #TASK 1c
    print(f"----------------------Task 1c----------------------")
    bearing_1 = UTM_bearing(latitude, longitude, radius(latitude, azimuth(local_dx1, local_dy1)), 0, local_dx1, 0, local_dy1, azimuth(local_dx1, local_dy1))
    bearing_2 = UTM_bearing(latitude, longitude, radius(latitude, azimuth(local_dx2, local_dy2)), 0, local_dx2, 0, local_dy2, azimuth(local_dx2, local_dy2))
    print(f"ST46-TP342 bearing: {bearing_1[0]} grad, compared to values from table 2: {azimuth(7033941.628-7034487.402, 568890.318 - 571578.304)} grad. Difference: {bearing_1[0] - azimuth(7033941.628-7034487.402, 568890.318 - 571578.304)} grad.")
    print(f"ST46-MOHOLT bearing: {bearing_2[0]} grad, compared to values from table 2: {azimuth(7031952.892-7034487.402, 571469.041 - 571578.304)} grad. Difference: {bearing_2[0] - azimuth(7031952.892-7034487.402, 571469.041 - 571578.304)} grad.")

    #TASK 1d
    print(f"----------------------Task 1d----------------------")
    NN2000_diff_1 = NN2000_height_diff(zenith(local_dx1, local_dy1, local_dz1), local_dx1, local_dy1, local_dz1, N_st46, N_tp342, radius(latitude, azimuth(local_dx1, local_dy1)))
    NN2000_diff_2 = NN2000_height_diff(zenith(local_dx2, local_dy2, local_dz2), local_dx2, local_dy2, local_dz2, N_st46, N_moholt, radius(latitude, azimuth(local_dx2, local_dy2)))
    print(f"ST46-TP342 NN2000 height difference: {NN2000_diff_1[0]} m, compared to values from table 2: {5.194-112.554} m. Difference: {NN2000_diff_1[0] - (5.194-112.554)}")
    print(f"ST46-MOHOLT NN2000 height difference: {NN2000_diff_2[0]} m, compared to values from table 2: {131.404-112.554} m. Difference: {NN2000_diff_2[0] - (131.404-112.554)}")
    #print(transformation_matrix(latitude, longitude, dx2, dy2, dz2))
    #print(distance_UTM((horizontal_distane(local_dx1, local_dy1)), , )
    #print(azimuth(local_dx1, local_dy1))
    #print(azimuth(local_dx2, local_dy2))


    #transform(transformation_matrix(test1, test2), dx1, dy1, dz1)

    print(f"----------------------Task 2a----------------------")
    covariance_1 = covariance_matrix(std_st46_tp342, K_st46_tp342)
    covariance_2 = covariance_matrix(std_st46_moholt, K_std46_moholt)
    print(f"Variance/Covariance matrix of st46-tp342 baseline: \n {covariance_1}")
    print(f"Variance/Covariance matrix of st46-moholt baseline: \n {covariance_2}")

    print(f"----------------------Task 2b----------------------")
    covariance_1_local = error_prop_law(transformation_matrix(latitude, longitude), covariance_1)
    covariance_2_local = error_prop_law(transformation_matrix(latitude, longitude), covariance_2)
    print(F_matrix(303.5776346, 9270.1001))
    covariance_1_obs = error_prop_law(F_matrix(azimuth(local_dx1, local_dy1), horizontal_distane(local_dx1, local_dy1)), covariance_1_local)
    covariance_2_obs = error_prop_law(F_matrix(azimuth(local_dx2, local_dy2), horizontal_distane(local_dx2, local_dy2)), covariance_2_local)
    print(np.sqrt(covariance_1_obs[0][0]), np.sqrt(covariance_1_obs[1][1]), np.sqrt(covariance_1_obs[2][2]))
    print(np.sqrt(covariance_2_obs[0][0]), np.sqrt(covariance_2_obs[1][1]), np.sqrt(covariance_2_obs[2][2]))


if __name__ == "__main__":
    main()