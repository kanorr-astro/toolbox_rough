#Simple Hohmann Transfer DV Calculator

#module imports
import math as m

def hohmann_dv(rad1, rad2, mu=1) ->float:
    vcs1 = m.sqrt(mu/rad1)
    print("VSC1: ", vcs1)
    vcs2 = m.sqrt(mu/rad2)
    print("VSC2:", vcs2)
    SME_transit = -mu / (rad1 + rad2)
    print("SME transit: ", SME_transit)
    vt1 = m.sqrt(2*((mu/rad1) + SME_transit))
    print("VT1: ", vt1)
    vt2 = m.sqrt(2*((mu/rad2) + SME_transit))
    print("VT2: ", vt2)
    dv1 = vt1 - vcs1
    print("DV1: ", dv1)
    dv2 = vcs2 - vt2
    print("DV2: ", dv2)
    return abs(dv1 + dv2)

if __name__ == "__main__":
    dv = hohmann_dv(3, 2)
    print("DeltaV Required for Transfer (DU/TU): ",dv)
