import numpy as np

t = 129600
we = 7.2921151467 * 10**(-5)
GM = 3.986005 * 10**14

def tk(toe):
    tk = t-toe
    if tk > 302400:
        tk -= 604800
    elif tk < -302400:
        tk += 604800
    return tk

def lambdak(lambda0, omega, toe):
    return lambda0 + (omega-we)*tk(toe) - we*toe

def Mk(a, deltan, M0, toe):
    return M0 + (np.sqrt(GM/a**6) + deltan) * tk(toe)

def EK(e, a, deltan, M0, toe, Ek = 1):
    Ek = 1
    for i in range(100):
        Ek = Mk(a, deltan, M0, toe) + e*np.sin(Ek)
    return Ek

def fk(e, a, deltan, M0, toe):
    return 2 * np.arctan(np.sqrt((1+e)/(1-e))*np.tan(EK(e, a, deltan, M0, toe)/2))

def ik(i0, i, cic, cis, w, toe, e, a, deltan, M0):
    return i0 + i*tk(toe) + cic*np.cos(2*(w+fk(e, a, deltan, M0, toe))) + cis*np.sin(2*(w+fk(e, a, deltan, M0, toe)))

def uk(w, e, a, deltan, M0, toe, cuc, cus):
    return w + fk(e, a, deltan, M0, toe) + cuc*np.cos(2*(w+fk(e, a, deltan, M0, toe))) + cus*np.sin(2*(w+fk(e, a, deltan, M0, toe)))

def rk(a, e, deltan, M0, toe, w, crc, crs):
    return (a**2)*(1-e*np.cos(EK(e, a, deltan, M0, toe))) + crc*np.cos(2*(w+fk(e,a,deltan,M0,toe))) + crs*np.sin(2*(w+fk(e,a,deltan,M0, toe)))

def R3(theta):
    return np.array([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])

def R1(theta):
    return np.array([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]])

toe = 1.295840*10**5
sqrta = 5.153681*10**3
e = 5.747278*10**-3
M0 = -2.941505
w = -1.770838
i0 = 9.332837*10**-1
lambda0 = 2.123898
deltan = 5.243075*10**-9
i = -6.853856*10**-10
omega = -8.116052*10**-9
cuc = -1.184642*10**-6
cus = 7.672235*10**-6
crc = 2.146562*10**2
crs = -2.140625*10**1
cic = 2.980232*10**-8
cis = -1.117587*10**-8

#print(np.array([[rk(sqrta, e, deltan, M0, toe, w, crc, crs)], [0], [0]]))

xk, yk, zk = R3(-lambdak(lambda0, omega, toe))@R1(-ik(i0, i, cic, cis, w, toe, e, sqrta, deltan, M0))@R3(-uk(w, e, sqrta, deltan, M0, toe, cuc, cus))@np.array([[rk(sqrta, e, deltan, M0, toe, w, crc, crs)], [0], [0]])
print(xk, yk, zk)

xk, yk, zk = R3(-lambdak(lambda0, 0, toe))@R1(-ik(i0, 0, 0, 0, w, toe, e, sqrta, 0, M0))@R3(-uk(w, e, sqrta, 0, M0, toe, 0, 0))@np.array([[rk(sqrta, e, 0, M0, toe, w, 0, 0)], [0], [0]])
print(xk, yk, zk)


#Task 3
a = 6378137
b = 6356752.3141

def N(a, b, lat):
    return a**2 / (np.sqrt(a**2*np.cos(lat)**2+b**2*np.sin(lat)**2))

lat = np.deg2rad(47.1)
long = np.deg2rad(15.5)
h = 400
print(np.array([(N(a,b,lat)+h)*np.cos(lat)*np.cos(long), (N(a,b,lat)+h)*np.cos(lat)*np.sin(long), ((b**2/a**2)*N(a,b,lat)+h)*np.sin(lat)]))