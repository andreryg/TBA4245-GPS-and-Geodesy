import numpy as np


class GPS(object):
    t = 129600
    we = 7.2921151467 * 10**-5
    GM = 3.986005 * 10**14
    
    def __init__(self, toe, sqrta, e, M0, w, i0, lambda0, deltan, i, omega, cuc, cus, crc, crs, cic, cis):
        self.toe = toe
        self.sqrta = sqrta
        self.e = e
        self.M0 = M0
        self.w = w
        self.i0 = i0
        self.lambda0 = lambda0
        self.deltan = deltan
        self.i = i
        self.omega = omega
        self.cuc = cuc
        self.cus = cus
        self.crc = crc
        self.crs = crs
        self.cic = cic
        self.cis = cis

    def disregard_corrections(self):
        self.deltan, self.i, self.omega, self.cuc, self.cus, self.crc, self.crs, self.cic, self.cis = [eval(i) for i in ["0"]*9]
        return True

    def _tk(self):
        tk = self.t-self.toe
        if tk > 302400:
            tk -= 604800
        elif tk < -302400:
            tk += 604800
        return tk

    def _lambdak(self):
        return self.lambda0 + (self.omega-self.we)*self._tk() - self.we*self.toe

    def _Mk(self):
        return self.M0 + (np.sqrt(self.GM/self.sqrta**6) + self.deltan) * self._tk()

    def _EK(self):
        Ek = self._Mk()
        for i in range(100):
            Ek = self._Mk() + self.e*np.sin(Ek)
        return Ek

    def _fk(self):
        return 2 * np.arctan(np.sqrt((1+self.e)/(1-self.e))*np.tan(self._EK()/2))

    def _ik(self):
        return self.i0 + self.i*self._tk() + self.cic*np.cos(2*(self.w+self._fk())) + self.cis*np.sin(2*(self.w+self._fk()))

    def _uk(self):
        return self.w + self._fk() + self.cuc*np.cos(2*(self.w+self._fk())) + self.cus*np.sin(2*(self.w+self._fk()))

    def _rk(self):
        return (self.sqrta**2)*(1-self.e*np.cos(self._EK())) + self.crc*np.cos(2*(self.w+self._fk())) + self.crs*np.sin(2*(self.w + self._fk()))

    def _R3(self, theta):
        return np.array([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])

    def _R1(self, theta):
        return np.array([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]])
    
    def calculateXYZ(self):
        return self._R3(-self._lambdak())@self._R1(-self._ik())@self._R3(self._uk())@np.array([[self._rk()], [0], [0]])

def N(lat, a=6378137, b=6356752.3141):
    return a**2 / (np.sqrt(a**2*np.cos(lat)**2+b**2*np.sin(lat)**2))


def main():
    #Task 1
    SV06 = GPS(1.295840*10**5, 5.153681*10**3, 5.747278*10**-3, -2.941505, -1.770838, 9.332837*10**-1, 2.123898, 5.243075*10**-9, -6.853856*10**-10, -8.116052*10**-9, -1.184642*10**-6, 7.672235*10**-6, 2.146562*10**2, -2.140625*10**1, 2.980232*10**-8, -1.117587*10**-8)
    SV10 = GPS(1.296000*10**5, 5.153730*10**3, 7.258582*10**-3, 4.044839*10**-1, 4.344642*10**-1, 9.71311*10**-1, -2.006987, 4.442685*10**-9, 2.521533*10**-10, -8.495353*10**-9, 4.714354*10**-6, -1.825392*10**-7, 3.868750*10**2, 8.978125*10**1, 3.725290*10**-9, 8.940696*10**-8)

    print(SV06.calculateXYZ())
    print(SV10.calculateXYZ())

    #Task 2
    SV06.disregard_corrections()
    print(SV06.calculateXYZ())


    #Task 3
    lat = np.deg2rad(47.1)
    long = np.deg2rad(15.5)
    a, b, h = 6378137, 6356752.3141, 400
    print(np.array([(N(lat)+h)*np.cos(lat)*np.cos(long), (N(lat)+h)*np.cos(lat)*np.sin(long), ((b**2/a**2)*N(lat)+h)*np.sin(lat)]))


if __name__=="__main__":
    main()