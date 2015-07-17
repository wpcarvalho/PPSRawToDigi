#!/usr/bin/python2

# authors: Jakub Sawicki <jakub.kuba.sawicki@gmail.com>
#          Jan Kaspar <jan.kaspar@gmail.com>

import re
import sys
from math import pi, cos, sin

# ie. usage:
# > cmsRun selected_events_analysis_cfg.py <module> 2 4 [RPs ...] | utils/verbose_to_asy.py [RPs ...] | asy -
# if your lxplus doesn't come with asymptote installed you can compile the output
# file at your local machine or install asymptote locally
#
# to run the asy file you will need pad_layout.asy file by Jan Kaspar<jan.kaspar@gmail.com>

RPs      = {}
dets     = {}
hits     = {}
clusters = []
tracks   = []
ev_id    = 0

interestingRPs = [int(RPId) for RPId in sys.argv[1:]]

def main():
    state = "begin"
    for line in sys.stdin:
        if line == "* RPs:\n":
            state = "RPs"
            continue
        elif line == "* silicon detectors:\n":
            state = "dets"
            continue
        elif re.match(r'>> RPFastStationSimulation::produce.*', line) != None:
            state = "sim"

            ev_id = int(re.match(r'>> RPFastStationSimulation::produce > event (\d+)\s*', line).group(1))
            continue
        elif re.match(r'>> RPNonParallelTrackCandidateFinder::produce.*', line) != None:
            state = "1rp_nonparallel"
            continue
        elif re.match(r'>> RPStationMultiTrackFinderFitter::produce.*', line) != None:
            state = "reco"
            continue
        elif re.match(r'>> RPStationMultiTrackFinderFitterAnalyzer::analyze.*', line) != None:
            state = "analyze"
            continue


        if state == "RPs":
            readRP(line)
        elif state == "dets":
            readDet(line)
        elif state == "sim":
            readSim(line)
        elif state == "1rp_nonparallel":
            read1RPNonParallel(line)
        elif state == "reco":
            readReco(line)
        elif state == "analyze":
            readAnalyze(line)

    patchRPInfo()

    generateCandidates()

    printAsyOutput(auto = True)


def readRP(line):
    # ID |   x (mm)    |   y (mm)    |   z  (m)    |
    #  0 |      +0.000 |    +107.825 | -148784.000
    match = re.match(r'\s*(\d+) \|\s*([+-]?\d+\.\d+) \|\s*([+-]?\d+\.\d+) \|\s*([+-]?\d+\.\d+)\s*', line)
    if match != None:
        if not int(match.group(1)) in interestingRPs:
            return

        RPs[int(match.group(1))] = { 'x' : float(match.group(2)),
                                     'y' : float(match.group(3)),
                                     'z' : float(match.group(4)) }


def readDet(line):
    #DetId |            detector center           |  readout direction
    #      |   x (mm)   |   y (mm)   |   z  (m)   |    dx    |    dy
    #   0  |    -0.514  |    54.717  |  -148.7637 |    0.707 |    0.707

    match = re.match(r'\s*(\d+)\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*([+-]?\d+\.\d+)\s*', line)
    if match != None:
        ID = int(match.group(1))

        if not ID/10 in interestingRPs:
            return

        dets[ID] = { 'x'         : float(match.group(2)),
                     'y'         : float(match.group(3)),
                     'z'         : float(match.group(4))*1e3,
                     'dx'        : float(match.group(5)),
                     'dy'        : float(match.group(6)),
                     'direction' : 'u' if ID % 2 != 0 else 'v' }


def patchRPInfo():
    '''
    Calculates the average direction of the u and v
    of each RP.
    '''
    for rpID, rpInfo in RPs.iteritems():
        u_dx = 0.
        u_dy = 0.
        u_count = 0

        v_dx = 0.
        v_dy = 0.
        v_count = 0

        for i in range(10):
            det = dets[rpID*10+i]
            if det['direction'] == 'u':
                u_dx += det['dx']
                u_dy += det['dy']
                u_count += 1
            else:
                v_dx += det['dx']
                v_dy += det['dy']
                v_count += 1

        rpInfo['u_dx'] = u_dx/u_count
        rpInfo['u_dy'] = u_dy/u_count
        rpInfo['v_dx'] = v_dx/v_count
        rpInfo['v_dy'] = v_dy/v_count


def readSim(line):
    pass


def read1RPNonParallel(line):
    match = re.match(r'\s*RP\s(\d+)\s*', line)
    if match != None:
        rpID = int(match.group(1))
        read1RPNonParallel.rpID = rpID
        hits[rpID] = { "u" : [], "v" : [] }
        return

    match = re.match(r'\s*u recognition\s*', line)
    if (match != None):
        read1RPNonParallel.direction = "u"
        return

    match = re.match(r'\s*v recognition\s*', line)
    if (match != None):
        read1RPNonParallel.direction = "v"
        return

    match = re.match(r'\s*a = ([-+]?\d+\.\d+), b = ([-+]?\d+\.\d+), w = ([-+]?\d+\.\d+)\s*', line)
    if (match != None):
        rpID = read1RPNonParallel.rpID
        direction = read1RPNonParallel.direction
        hits[rpID][direction].append({
            'a' : float(match.group(1)),
            'b' : float(match.group(2)),
            'w' : float(match.group(3)) })

read1RPNonParallel.rpId = 0
read1RPNonParallel.direction = ""


def readReco(line):
    # * selected cluster #0: ax = +10.0 urad, ay = +1.2 urad, bx = -0.297 mm, by = +4.519 mm, z0 = 212.698 m
    match = re.match(r'\* selected cluster #(\d+): ax = ([-+]?\d+\.\d+) urad, ay = ([-+]?\d+\.\d+) urad, bx = ([-+]?\d+\.\d+) mm, by = ([-+]?\d+\.\d+) mm, z0 = ([-+]?\d+\.\d+) m\s*', line)
    if match != None:
        clusters.append({
            'ax'     : float(match.group(2))*1e-6,
            'ay'     : float(match.group(3))*1e-6,
            'bx'     : float(match.group(4)),
            'by'     : float(match.group(5)),
            'z0'     : float(match.group(6))*1e3,
            'chosen' : True })
        return

    # get info about track association with the hits
    match = re.match(r'cluster #(\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+)\s*', line)
    if match != None:
        clusters[int(match.group(1))]['hits'] = {
                100 : { 'u' : int(match.group(2)), 'v' : int(match.group(3)) },
                104 : { 'u' : int(match.group(4)), 'v' : int(match.group(5)) },
                120 : { 'u' : int(match.group(6)), 'v' : int(match.group(7)) } }
        return

    # has the track been chosen
    match = re.match(r'(\d+) (true|false)\s*', line)
    if match != None:
        clusters[int(match.group(1))]['chosen'] = True if match.group(2) == 'true' else False


def readAnalyze(line):
    match = re.match(r'\s*ax=([-+]?\d+\.\d+) urad, ay=([-+]?\d+\.\d+) urad, bx=([-+]?\d+\.\d+) mm, by=([-+]?\d+\.\d+) mm, z0=([-+]?\d+\.\d+) mm, reconstructable=1\s*', line)
    if match != None:
        tracks.append({
            'ax' : float(match.group(1))*1e-6,
            'ay' : float(match.group(2))*1e-6,
            'bx' : float(match.group(3)),
            'by' : float(match.group(4)),
            'z0' : float(match.group(5)) })


def generateCandidates():
    '''
    Generate clusters used for Hough transform.
    Combines hits in u and v directions to get possible hit points.
    '''
    for RPID, RPInfo in RPs.iteritems():
        z = RPInfo['z']
        us = []
        vs = []

        u_dx = RPInfo['u_dx']
        u_dy = RPInfo['u_dy']
        v_dx = RPInfo['v_dx']
        v_dy = RPInfo['v_dy']

        # get the simulated potential clusters
        for track in tracks:
            z0 = track['z0']
            x = track['ax'] * (z - z0) + track['bx']
            y = track['ay'] * (z - z0) + track['by']

            u = u_dx * x + u_dy * y
            v = v_dx * x + v_dy * y

            us.append(u)
            vs.append(v)

        cand = []

        for u in us:
            for v in vs:
                cand.append((v_dy * u - u_dy * v, -v_dx * u + u_dx * v))

        RPInfo['candidates_real'] = cand

        # get the actual reconstructed clusters used in further reconstruction
        cand = []
        hitz = hits[RPID]

        for u in hitz['u']:
            for v in hitz['v']:
                cand.append((v_dy * u['b'] - u_dy * v['b'], -v_dx * u['b'] + u_dx * v['b']))

        RPInfo['candidates_reco'] = cand


def printAsyOutput(auto = False, center = (0., 4.5), width = 2., z_min = 202e3, z_max = 222e3):
    print "import pad_layout;"

    if auto:
        x_max = -1e100
        y_max = -1e100
        x_min = +1e100
        y_min = +1e100

        for rpID, rpInfo in RPs.iteritems():
            for x, y in rpInfo['candidates_real'] + rpInfo['candidates_reco']:
                if x > x_max:
                    x_max = x
                if y > y_max:
                    y_max = y
                if x < x_min:
                    x_min = x
                if y < y_min:
                    y_min = y

        width = max(x_max - x_min, y_max - y_min) * 1.3
        center = ((x_max + x_min) / 2, (y_max + y_min) / 2 + width / 1.3 * 0.1)

    top    = center[1] + width / 2
    bottom = center[1] - width / 2
    right  = center[0] + width / 2
    left   = center[0] - width / 2

    # XY projection of active strips and reconstructed tracks' hits
    for rpID in interestingRPs:
        rpInfo = RPs[rpID]
        z = rpInfo['z']

        print "NewPad(\"$x \ung{mm}$\", \"$y \ung{mm}$\");"

        # strips
        for direction in ['u', 'v']:
            hitz = hits[rpID][direction]
            for m in hitz:
                x_min = left
                x_max = right

                if rpInfo[direction+"_dy"] == 0:
                    continue

                # perpendicular to the readout direction lies the strip
                a = - rpInfo[direction+"_dx"] / rpInfo[direction+"_dy"]
                b = m['b'] / rpInfo[direction+"_dy"]

                y_min = a * x_min + b
                y_max = a * x_max + b
                print "draw(({0},{1})--({2},{3}),{4});".format(x_min,y_min,x_max,y_max,"solid+lightgrey+3pt")

        # actual hits
        i = 1
        for track in tracks:
            x = track['ax'] * (z - track['z0']) + track['bx']
            y = track['ay'] * (z - track['z0']) + track['by']

            marker = "mCi+1pt+StdPen({0})".format(i)
            print "draw(({0},{1}),{2});".format(x, y, marker)
            i += 1

        # reconstructed hits
        i = 1
        for track in clusters:
            x = track['ax'] * (z - track['z0']) + track['bx']
            y = track['ay'] * (z - track['z0']) + track['by']

            marker = "mCi+3pt+StdPen({0})+false".format(i)
            print "draw(({0},{1}),{2});".format(x, y, marker)
            i += 1

        print "limits(({0},{1}),({2},{3}),Crop);".format(left, bottom, right, top)
        print "AttachLegend(\"RP{0}, z = {1}\");".format(rpID, z)

    print "NewRow();"

    # XY projection of hits and clusters in each RP plane
    for rpID in interestingRPs:
        rpInfo = RPs[rpID]
        z = rpInfo['z']

        print "NewPad(\"$x \ung{mm}$\", \"$y \ung{mm}$\");"

        # reconstructed hits
        i = 1
        for track in clusters:
            x = track['ax'] * (z - track['z0']) + track['bx']
            y = track['ay'] * (z - track['z0']) + track['by']

            marker = "mCi+3pt+false+StdPen({0})".format(i)
            print "draw(({0},{1}),{2});".format(x, y, marker)
            i += 1

        # clusters
        i = 1
        for x, y in rpInfo['candidates_reco']:
            print "draw(({0},{1}),mCr+3pt+black);".format(x, y)
            print "label(\"{0}\",({1},{2}),E);".format(str(i), x, y)
            i += 1

        # actual hits
        i = 1
        for track in tracks:
            x = track['ax'] * (z - track['z0']) + track['bx']
            y = track['ay'] * (z - track['z0']) + track['by']

            marker = "mCi+1pt+StdPen({0})".format(i)
            print "draw(({0},{1}),{2});".format(x, y, marker)
            i += 1

        print "limits(({0},{1}),({2},{3}),Crop);".format(left, bottom, right, top)
        print "AttachLegend(\"RP{0}, z = {1}\");".format(rpID, z)

    print "NewRow();"

    # ZX and ZY projections of tracks and clusters parallel to the beam direction
    for angle in [0, pi/2]:
        print "NewPad(\"z [m]\",\"projection in ({0:.2f}, {1:.2f}) [mm]\");".format(cos(angle), sin(angle))
        dx = cos(angle)
        dy = sin(angle)

        # clusters
        for rpID, rpInfo in RPs.iteritems():
            z = rpInfo['z']*1e-3

            i = 1
            for x, y in rpInfo['candidates_reco']:
                p = x * dx + y * dy
                print "draw(({0},{1}),{2});".format(z,p,"mCr+black+3pt")
                print "draw(({0}, {1}-0.033)--({0}, {1}+0.033), black);".format(z, p)
                print "label(\"{0}\",({1},{2}),E);".format(str(i), z, p)
                i += 1

        # tracks
        for tracks_, pen in [(tracks, ""),(clusters, "+dashed")]:
            i = 1
            for track in tracks_:
                x_min = track['ax'] * (z_min - track['z0']) + track['bx']
                y_min = track['ay'] * (z_min - track['z0']) + track['by']

                x_max = track['ax'] * (z_max - track['z0']) + track['bx']
                y_max = track['ay'] * (z_max - track['z0']) + track['by']

                p_min = x_min * dx + y_min * dy
                p_max = x_max * dx + y_max * dy

                try:
                    p = pen if track['chosen'] else "+Dotted"
                except KeyError:
                    p = pen
                print "draw(({0},{1})--({2},{3}),StdPen({4}){5});".format(z_min*1e-3, p_min, z_max*1e-3, p_max, i, p)
                i += 1

    # put some legend for non magicians to see clearer
    print "NewPad(false);"
    for items in [
            ["\"<true tracks (each track different colour):\""],
            ["\"hit\"", "mCi+1pt+StdPen(1)"],
            ["\"track\"", "StdPen(1)"],

            ["\"<hit candidates:\""],
            ["\"active strip\"", "solid+lightgrey+3pt"],
            ["\"U+V strip intersections (numbered)\"", "mCr+black+3pt"],

            ["\"<track candidatess (each candidate different colour):\""],
            ["\"any TC (in views along beam)\"", "mCi+3pt+StdPen(1)+false"],
            ["\"accepted TC\"", "StdPen(1)+dashed"],
            ["\"rejected TC\"", "StdPen(1)+Dotted"]]:
        print "AddToLegend({0});".format(",".join(items))
    print "AttachLegend(legAlig=(0,0),picAlig=(0,0));"


if __name__ == '__main__':
    main()
