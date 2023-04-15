#Simple Hohmann Transfer DV Calculator

#module imports
import math as m

def hohmann_dv(rad1, rad2, mu=1) ->float:
    vcs1 = m.sqrt(mu/rad1)
    vcs2 = m.sqrt(mu/rad2)
    SME_transit = -mu / (rad1 + rad2)
    v1 = m.sqrt(2*((mu/rad1) + SME_transit))
    v2 = m.sqrt(2*((mu/rad2) + SME_transit))
    dv1 = v1 - vcs1
    dv2 = v2 - vcs2
    return dv1 + dv2

if __name__ == "__main__":
    dv = hohmann_dv(2,3)
    print(dv)
