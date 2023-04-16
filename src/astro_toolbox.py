#Astrodynamics Toolbox

#module imports
import math as m
import numpy as np

#General Tools

def magnitude(vector):
    #vector magnitude
    return m.sqrt(sum(pow(element, 2) for element in vector))

def spec_mech_energy(r,v,mu):
    #calculate specific mechanical energy
    SME = (magnitude(v)**2 / 2) - (mu / magnitude(r))
    return SME   

def radius(p, e, nu):
    #calculate orbital radius for a given p,e, and true anomaly(nu)
    return p / (1 + e * m.cos(nu))

def axis_calcs(SME, h, mu=1):
    #calculate semi-major axis and semi-latus rectum
    if SME >= 0:
        a = 'inf'
        p = magnitude(h)**2 / mu
    else:
        a = -mu / (2*SME)
        p = magnitude(h)**2 / mu
    return a,p

def conic_type(h,e) ->str:
    #Classify conic type from angular momentum (h) and eccentricity (e)
    if magnitude(h) > 0:
        if e == 0:
            conic = "Circular"
        elif (e > 0) and (e < 1):
            conic = "Elliptical"
        elif e == 1:
            conic = "Parabolic"
        elif e > 1:
            conic = "Hyperbolic"
    else:
        conic = "Degenerate"
    return conic
    

#Orbital Determination and Coordinate Changes

def rv_orbparams(ri, rj, rk, vi, vj, vk, disp_params=0, mu=1) ->list:
    #Calculates orbital params for given r and v vectors
    r = np.array([ri, rj, rk])
    v = np.array([vi, vj, vk])

    #Determine orbital vectors (h, n, e)
    h_vect = np.cross(r,v)
    n_vect = np.cross([0,0,1],h_vect)
    e_vect = (1/mu)*(((magnitude(v)**2 - (mu/magnitude(r)))*r) - (np.dot(r,v)*v))

    #Determine orbital params
    SME = spec_mech_energy(r,v,mu)
    a,p = axis_calcs(SME, h_vect, mu)
    e = m.sqrt(1 + ((2*SME*magnitude(h_vect)**2)/(mu**2)))
    conic = conic_type(h_vect,e)
    i = m.acos(h_vect[2] / magnitude(h_vect))
    ascend = m.acos(n_vect[0] / magnitude(n_vect))
    peri = m.acos(np.dot(n_vect,e_vect) / (magnitude(n_vect)*magnitude(e_vect)))

    #display parameters (if selectred)
    if disp_params == 1:
        print("Orbit Conic Type: ", conic)
        print("Semi-Major Axis(a): ", a)
        print("Semi-Latus Rectum(p): ", p)
        print("Eccentricity(e): ", e)
        print("Inclination(i): ", m.degrees(i))
        print("Longitude of Ascending Node(Omega): ", m.degrees(ascend))
        print("Argument of Periapsis: ", m.degrees(peri))    
    return (a, p, e, i, ascend, peri, SME, conic)

def param_to_pqw(p, e, points=100) ->list:
    #Generate anomaly angles to evaluate (based on conic type)
    nu_raw = np.linspace(0,2*m.pi, points)
    nu_mod = []
    if e >= 1.25:
        for nu in nu_raw:
            if (nu > (m.pi / 2)) and (nu < (3*m.pi / 2)):
                None
            else:
                nu_mod.append(nu)
    elif (e > 1) and (e < 1.25):
        for nu in nu_raw:
            if (nu > (3*m.pi / 4)) and (nu < (5*m.pi / 4)):
                None
            else:
                nu_mod.append(nu)
    else:
        nu_mod = nu_raw
    
    #Using nu angles generate list of points in PQW orientation
    PQW = []
    for nu in nu_mod:
        PQW.append(np.array([radius(p,e,nu)*m.cos(nu), radius(p,e,nu)*m.sin(nu), 0]))
    return PQW

def pqw_to_ijk(pqw, i, ascend, peri) ->list:
    #define rotation matrix
    R_raw = [m.cos(ascend)*m.cos(peri)-m.sin(ascend)*m.sin(peri)*m.cos(i),
         -m.cos(ascend)*m.sin(peri)-m.sin(ascend)*m.cos(peri)*m.cos(i),
         m.sin(ascend)*m.sin(i),
         m.sin(ascend)*m.cos(peri)+m.cos(ascend)*m.sin(peri)*m.cos(i),
         -m.sin(ascend)*m.sin(peri)+m.cos(ascend)*m.cos(peri)*m.cos(i),
         -m.cos(ascend)*m.sin(i),
         m.sin(peri)*m.sin(i),
         m.cos(peri)*m.sin(i),
         m.cos(i)]
    R = np.array(R_raw).reshape(3,3)

    #Transform from PQW to IJK
    IJK = []
    for point in pqw:
        IJK.append(np.dot(R,point))
    return IJK


#Orbital Manuever Calculation Tools

def hohmann_dv(rad1, rad2, dv_disp=0, mu=1) ->float:
    #Calculates deltaV required for Hohmann Transfer between circular orbits of radius1 and radius2
    vcs1 = m.sqrt(mu/rad1)
    vcs2 = m.sqrt(mu/rad2)
    SME_transit = -mu / (rad1 + rad2)
    vt1 = m.sqrt(2*((mu/rad1) + SME_transit))
    vt2 = m.sqrt(2*((mu/rad2) + SME_transit))
    dv1 = vt1 - vcs1
    dv2 = vcs2 - vt2
    if dv_disp == 1:
        print("VCS1: ", vcs1)
        print("VCS2: ", vcs2)
        print("SME Tranist Orbit: ", SME_transit)
        print("VT1: ", vt1)
        print("VT2: ", vt2)
        print("DV1: ", dv1)
        print("DV2: ", dv2)
    return abs(dv1 + dv2)