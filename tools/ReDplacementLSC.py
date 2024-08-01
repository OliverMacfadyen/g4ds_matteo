#!/usr/bin/python

import ast
import getopt
import math
import sys
import numpy

cmPrint = False
#cmPrint = True

# units
m   = 1.
cm  = m / 100.
s   = 1.
ns  = 1E-9 * s
MeV = 1.
GeV = 1000. * MeV
keV = MeV / 1000.0

# constants
nMass = 939.565 * MeV
Ar40Mass = 37.2244 * GeV
c = 299792458. * m / s

X = -1.38496  * m
Y =  0.313385 * m
Z =  0.483399 * m
T = 8.1 * MeV
TH_N =  18.8   # angle vs. xy plane
PH_N = -12.75  # angle vs. x-axis
TH_S =  37.6   # scattering angle

TPC_POS = numpy.array( [ 0, 0, 0 ] )
D = 80 * cm    # distance TPC - LSci
N = 9          # number of LSci
PH_MIN =  90   # angle of the first LSC in plane perpendicular to neutron (90 points up)
PH_MAX = 180   # angle of the last LSC (180 points right)

gun_x = X
gun_y = Y
gun_z = Z
gun_t = T
th_n = TH_N
ph_n = PH_N
th_s = TH_S
tpc_pos = TPC_POS
d = D
n = N
ph_min = PH_MIN
ph_max = PH_MAX
verbose = False

def usage():
    #      12345678901234567890123456789012345678901234567890123456789012345678901234567890
    print 'Usage: ' + sys.argv[0] + ' [OPTION]...'
    print ''
    print 'Script to create a text file with the positions of the LSci, to be read by g4ds.'
    print ''
    print '  -x, --gun_x=VALUE    gun x position, in cm (default', X, 'm)'
    print '  -y, --gun_y=VALUE    gun y position, in cm (default', Y, 'm)'
    print '  -z, --gun_z=VALUE    gun z position, in cm (default', Z, 'm)'
    print '  -e, --gun_e=VALUE    gun kinetic energy, in keV (default', T, 'MeV)'
    print '  -t, --theta_n=NUM    neutron incident angle in xz plane (default %.2f),' % TH_N
    print '                       positive angle goes below horizontal (xy) plane'
    print '  -p, --phi_n=NUM      neutron incident angle vs. xz plane (default %.2f),' % PH_N
    print '  -s, --theta_s=NUM    neutron scattering angle (default %.2f)' % TH_S
    print '  -d, --distance=NUM   LSci distance to TPC center, in cm (default', D, 'cm)'
    print '  -#, --num=NUM        number of LSci (default', N, ')'
    print '  -m, --phi_min=VALUE  minimum azimutal angle (default', PH_MIN, ')'
    print '  -M, --phi_max=VALUE  maximum azimutal angle (default', PH_MAX, ')'
    print '  -v, --verbose        print extra info, to debug the LSci placement'
    print '  -h, --help           display this help and exit'
    print ''
    print 'Report bugs to <Michael.Kuss@pi.infn.it>.'

try:
    opts, args = getopt.getopt(sys.argv[1:], 'hx:y:z:e:t:p:s:d:#:m:M:v', [ 'help', 'theta_n=', 'phi_n=', 'theta_s=', 'distance=', 'num=', 'phi_min=', 'phi_max=', 'verbose' ])
except getopt.GetoptError, err:
    print err
    usage()
    sys.exit(2)
for o, a in opts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-x', '--gun_x'):
        gun_x = float(a)
    elif o in ('-y', '--gun_y'):
        gun_y = float(a)
    elif o in ('-z', '--gun_z'):
        gun_z = float(a)
    elif o in ('-e', '--gun_e'):
        gun_t = float(a)
    elif o in ('-t', '--theta_n'):
        th_n = float(a)
    elif o in ('-p', '--phi_n'):
        ph_n = float(a)
    elif o in ('-s', '--theta_s'):
        th_s = float(a)
    elif o in ('-d', '--distance'):
        d = float(a)
    elif o in ('-#', '--num'):
        n = int(a)
    elif o in ('-m', '--phi_min'):
        ph_min = float(a)
    elif o in ('-M', '--phi_max'):
        ph_max = float(a)
    elif o in ('-v', '--verbose'):
        verbose = True
    else:
        assert False, 'unhandled option'

gun_pos = [ gun_x, gun_y, gun_z ]
sTh_n = math.sin(math.radians(th_n))
cTh_n = math.cos(math.radians(th_n))
sPh_n = math.sin(math.radians(ph_n))
cPh_n = math.cos(math.radians(ph_n))
sTh_s = math.sin(math.radians(th_s))
cTh_s = math.cos(math.radians(th_s))

# rotation around z axis by neutron phi
rotnz = numpy.array( [ [ cPh_n, -sPh_n, 0 ],
                       [ sPh_n,  cPh_n, 0 ],
                       [     0,      0, 1 ] ] )
# rotation around y axis by neutron theta
rotny = numpy.array( [ [  cTh_n, 0, sTh_n ],
                       [      0, 1,      0 ],
                       [ -sTh_n, 0, cTh_n ] ] )
rotn = numpy.dot(rotnz, rotny)

if n > 1:
    dPhi = ( math.radians(ph_max) - math.radians(ph_min) ) / (n-1)
else:
    dPhi = 0

# rotation around z axis by scattering angle
rotzs = numpy.array( [ [ cTh_s, -sTh_s, 0 ],
                         [ sTh_s,  cTh_s, 0 ],
                         [     0,      0, 1 ] ] )

LSciPos = []

for i in xrange(n):
    phi = math.radians(ph_min) + dPhi * i;
    sPhi = math.sin(phi)
    cPhi = math.cos(phi)
    # rotation around x axis by phi, create individual position of LSci on ring
    rotPhi = numpy.array( [ [ 1,    0,     0 ],
                            [ 0, cPhi, -sPhi ],
                            [ 0, sPhi,  cPhi ] ] )
    # rotation around y axis by -theta_n (minus because xz is left handed)
    rotN = numpy.array( [ [  cTh_n, 0, sTh_n ],
                          [      0, 1,      0 ],
                          [ -sTh_n, 0, cTh_n ] ] )
    rot = numpy.dot(rotn, numpy.dot(rotPhi, rotzs))
    pos = numpy.dot(rot, [ d, 0, 0 ])
    LSciPos.append(pos)


nDir = tpc_pos - gun_pos
gunTPCdist = numpy.linalg.norm(nDir)
nDirUnit = nDir / gunTPCdist
nAng = numpy.arccos( numpy.dot( nDirUnit, [1,0,0] ) )
nMomNorm = math.sqrt( 2.0 * nMass * gun_t )
nVelNorm = nMomNorm / nMass
nMomVec  = nMomNorm * nDirUnit
nVelVec  = nVelNorm * nDirUnit
gunTPCToF = gunTPCdist / nVelNorm / c # neutron ToF gun - TPC
if verbose:
    print '# gun:'
    print '#   position %.3f %.3f %.3f m ==> distance %.3f m  angle %.1f deg' % ( gun_pos[0]/m, gun_pos[1]/m, gun_pos[2]/m, gunTPCdist/m, math.degrees(nAng) )
    print '#   T %.3f MeV ==> p %.6f MeV/c v/c %.3f [ %.3f %.3f %.3f ] c  ToF %.1f ns' % ( gun_t/MeV, nMomNorm/MeV, nVelNorm, nVelVec[0], nVelVec[1], nVelVec[2], gunTPCToF/ns )
    print '#'

# cm
cmMass = nMass + Ar40Mass
cmVelVec = (nMass * nVelVec + Ar40Mass * 0 ) / cmMass
cmVelNorm = numpy.linalg.norm(cmVelVec)
cmE = cmMass + gun_t
cmT = 0.5 * cmMass * cmVelNorm * cmVelNorm
cmEx = cmE - cmMass - cmT
nVelNormCM = nVelNorm - cmVelNorm
Ar40VelNormCM = cmVelNorm 
if verbose:
    print '# cm: speed %.5f c [ %.5f %.5f %.5f ] c  T %.4f MeV  Ex %.4f MeV' % ( cmVelNorm, cmVelVec[0], cmVelVec[1], cmVelVec[2], cmT/MeV, cmEx/MeV )
    print '#'

if verbose:
    print '# LSci\'s:'
for pos in LSciPos:
    n2Dir = pos - tpc_pos
    TPCLSciDist = numpy.linalg.norm(n2Dir)
    n2DirUnit = n2Dir / TPCLSciDist
    n2Ang = numpy.arccos( numpy.dot( nDirUnit, n2DirUnit ) )
    print '%6.2f cm  %6.2f cm  %6.2f cm' % ( pos[0] / cm, pos[1] / cm, pos[2] / cm )
    print '#   distance %6.3f m scattAngle %.3f deg' % ( TPCLSciDist/m, math.degrees(n2Ang) )
    # which angle does the neutron have in cm (wrt to nDir) to have the desired lab scattering angle?
    # this all is now in the plane spanned by n_dir and n2_dir
    #
    #                   _-/
    #                 _-A/
    #               _-  /
    #             _-   /
    # v_n_lab=c _-    / b=v_n_cm      x/sin(X) = const
    #         _-     /
    #       _-      /
    #     _- B   C /
    #    ---------/
    #        a = v_cm

    lawOfSines = nVelNormCM / math.sin(n2Ang)
    gamma = math.asin( cmVelNorm / lawOfSines )
    n2Angcm = gamma + n2Ang  # complementary of 180 - gamma - n2Ang
    Ar40Angcm = math.pi + n2Angcm
    if verbose:
        print '#   cm:  ',
        print 'neutron: v/c %.6f angle %6.2f deg  ' % ( nVelNormCM,  math.degrees(n2Angcm) ),
        print 'Ar40: v/c %.5f  angle %6.2f deg' % ( Ar40VelNormCM, math.degrees(Ar40Angcm) )
    # simple calculation:
    # 1) get the momentum of the neutron versus the LSci
    # 2) determine the momentum vector of the Ar40
    # 3) check if kinetic energy is conserved
    n2VelVec = lawOfSines * math.sin(n2Angcm)  # 180-ang, but for the sin it's the same
    n2VelNorm = numpy.linalg.norm(n2VelVec)
    n2MomVec = n2VelVec * n2DirUnit * nMass
    n2T = 0.5 * nMass * n2VelVec * n2VelVec
    TPCLSciToF = TPCLSciDist / n2VelNorm / c # neutron ToF TPC -> LSci
    print '#   lab: neutron: v/c %.4f T %.6f MeV  ToF %.1f ns' % ( n2VelVec, n2T/MeV, TPCLSciToF/ns )
    Ar40MomVec = nMomVec - n2MomVec
    Ar40MomNorm = numpy.linalg.norm(Ar40MomVec)
    Ar40MomUnit = Ar40MomVec / Ar40MomNorm
    Ar40VelVec = Ar40MomVec / Ar40Mass
    Ar40VelNorm = numpy.linalg.norm(Ar40MomVec) / Ar40Mass
    Ar40T = 0.5 * Ar40Mass * Ar40VelNorm * Ar40VelNorm
    Ar40cosTh = numpy.dot( [0,0,1], Ar40MomUnit )
    Ar40Th = math.degrees(math.acos(Ar40cosTh))
    print '#   lab: Ar40: v/c %6.4f (%7.4f,%7.4f,%7.4f)  th %3.0f (%5.2f)  T %6.3f keV' % ( Ar40VelNorm, Ar40VelVec[0], Ar40VelVec[1], Ar40VelVec[2], Ar40Th, Ar40cosTh, Ar40T/keV )
    ratioT = gun_t / ( n2T + Ar40T ) - 1.0
    print '#   debug: ratioT %g' % ( ratioT )
    if cmPrint:
        break
