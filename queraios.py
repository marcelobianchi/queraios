#! /usr/bin/python

################################################################################
#  QueRaios: an evaluation tool of crossing raypaths density in                #
#  Ambient Noise Tomography.                                                   #
#                                                                              #
#   Copyright (C) 2016  Marcelo B. Bianchi, Bruno Collaco                      #
#                                                                              #
#  This program is free software: you can redistribute it and/or modify        #
#  it under the terms of the GNU General Public License as published by        #
#  the Free Software Foundation, either version 3 of the License, or           #
#  (at your option) any later version.                                         #
#                                                                              #
#  This program is distributed in the hope that it will be useful,             #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of              #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
#  GNU General Public License for more details.                                #
#                                                                              #
#  You should have received a copy of the GNU General Public License           #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.       #
################################################################################

#
# Grid resolution computation by considering combinations of rays 2by2 rays
#

from __future__ import print_function
from __future__ import division

from scipy.io import netcdf_file as netcdf
from scipy.misc import comb
from scipy import version
import numpy
import glob
import sys
import os
import shutil
import resource

from optparse import OptionParser
from optparse import OptionGroup

def using(point=""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: User time=%s System time=%s Memory=%.2f mb''' % (
        point,
        usage[0],
        usage[1],
        (usage[2] * resource.getpagesize()) / 1000000.0)


def fitnumber(number, res, right=False):
    v = res * (number // res)
    if right and v < number:
        v += res
    return v


def parserays(filename):
    rays = []
    lines = open(filename)
    for line in lines:
        items = line.strip().split()
        items = list(map(float, items))
        rays.append((items[0], items[1], items[2], items[3]))
    return rays


def write(filename, data):
    nfile = netcdf(filename, 'a')
    nfile.variables['z'][::] = data
    nfile.close()
    return


def read(filename):
    nfile = netcdf(filename, 'r')
    data = nfile.variables['z'][::]
    data = data.copy()
    nfile.close()
    return numpy.nan_to_num(data)


def analyse(a, b, r):
    c = (a + b)
    c[c == 1] = 0  # Remove non overlaps
    c[c == 2] = 1  # Normalize to 1
    nca = numpy.sum(a)
    ncb = numpy.sum(b)
    ncj = numpy.sum(c)
    if ncj == 0:
        return r, 0
    norm = ((nca - ncj + 1) * (ncb - ncj + 1)) / (nca * ncb)
    r += c * norm
    return r, 1


def combinereal(this, gds):
    r = (gds[this]).copy() * 0.0

    cc = 0
    jj = 0
    for that in range(this + 1, len(gds)):
        r, j = analyse(gds[this], gds[that], r)
        jj += j
        cc += 1

    print("\rDone processing for Ray %05d = " % this, end="")
    return cc, jj, r


def combine(gds):
    tc = comb(len(gds), 2)  # Total combinations possible
    cc = 0                  # Current combination
    jj = 0                  # Total considered combination

    rr = gds[0].copy() * 0.0
    for this in range(len(gds)):
        c, j, r = combinereal(this, gds)
        rr += r
        jj += j
        cc += c
        print ('%09d [%03d %%]' % (jj, 100 * cc / tc), end="")

    print("\n Considered %d of %d rays crossing paths." % (jj, tc))
    rr = rr / jj
    rr[rr == 0] = numpy.nan
    return rr


def hitcount(gds):
    final = gds[0] * 0.0
    for this in gds:
        final += this
    final /= len(gds)
    return final


def buildgrids(rays, limits, sps, resolution, basename, dirname, recreate=False, length=False):
    R = "-R%f/%f/%f/%f" % limits
    I = "-I%f/%f"       % resolution
    i = 0
    for lon1, lat1, lon2, lat2 in rays:
        i += 1
        print("Computing grids %05d/%05d = %d %%" % (i, len(rays), 100 * i / len(rays)), end="\r")

        filename = os.path.join(dirname,"data-%05d-%s.grd" % (i, basename))
        if recreate is False and os.path.isfile(filename):
            continue

        if length:
            # Computing ray path in km inside box - needs to be checked
            cmd = "project -C%f/%f -E%f/%f -Q -G%f | awk '{print $1,$2,1}' | blockmean -Sz -C %s %s | awk -v sps=%f '{print $1,$2,$3*sps}' | xyz2grd %s %s -G%s" % (lon1, lat1, lon2, lat2, sps, R, I, sps, R, I, filename)
        else:
            cmd = "project -C%f/%f -E%f/%f -Q -G%f | awk '{print $1,$2,1}' | blockmean -C %s %s | xyz2grd %s %s -G%s" % (lon1, lat1, lon2, lat2, sps, R, I, R, I, filename)

        os.system(cmd)
    print()


def plotresults(datafile, rayfile, psfile, maxcolor, limits):
    xinc = "%.3f" % (maxcolor / 5)
    title = "Rays Cross Indicator (High is better)"
    R = "-R%f/%f/%f/%f" % limits

    os.system("gmtset PAPER_MEDIA a0")
    os.system("grd2cpt   %s -L0/%f -Z -M  > la.cpt" % (datafile, maxcolor))

    os.system("psbasemap %s -B15WsNE -JM11 -X3 -Yc -K > %s"                 % (R, psfile))
    os.system("grdview   %s -To -Cla.cpt -R -J -O -K >> %s"                 % (datafile, psfile))
    os.system("pscoast   -R -J -O -K -Na -W1p,black -A10000 >> %s"          % (psfile))
    os.system('psscale   -D6/-1/10/0.5h -O -K -Cla.cpt -B%s:."%s": -E >> %s' % (xinc, title, psfile))
    os.system("awk -v V='\n' '{print $1,$2,V$3,$4}' %s | sort -u | psxy -R -J -O -K -Sc0.2 -G200 -m >> %s" % (rayfile, psfile))
    # os.system("awk -v V='\n' '{print \">\"V$1,$2,V$3,$4}' %s | psxy -R -J -O -K -W0.5p -m >> %s" % (rayfile, psfile))

    os.system("psbasemap -R -B15WsNE -J -X14  -K -O >> %s"                  % (psfile))
    os.system("pscoast   -R -J -O -K -Na -W1p,black -A10000 >> %s"          % (psfile))
    os.system("awk -v V='\n' '{print \">\"V$1,$2,V$3,$4}' %s | psxy -R -J -O -K -W1p,gray -m >> %s"           % (rayfile, psfile))
    os.system("echo '' | psxy -R -J -O >> %s"                               % (psfile))

    os.system("ps2raster -P -A -Tf %s"                                      % (psfile))
    os.system("rm -f la.cpt %s"                                             % (psfile))


def checkSystem():
    tools = [
        "gmtset",
        "grd2cpt",
        "psbasemap",
        "grdview",
        "pscoast",
        "psscale",
        "awk",
        "sort",
        "psxy",
        "project",
        "blockmean",
        "xyz2grd"
    ]

    top, middle, smaller = version.version.split(".")
    if int(middle) < 15:
        print("** Scipy version should be >= 0.15.0, current version is %s **" % version.version)
        sys.exit(1)

    missing = []

    for t in tools:
        r = os.system("which %s > /dev/null 2>&1" % t)
        if r != 0:
            missing.append(t)
    return missing


if __name__ == "__main__":
    # Check that this system has all the necessary tools for using this code
    #
    missingtools = checkSystem()
    if missingtools:
        print("Necessary tool%s for running this code %s missing:\n  %s" %
              (("" if len(missingtools) == 1 else "s"), ("are" if len(missingtools) != 1 else "is"), missingtools))
        sys.exit(1)

    # Build command line arguments
    #
    parser = OptionParser()
    
    main = OptionGroup(parser, "Main Mandatory Options", "Those options are mandatory for the code to work.")
    opti = OptionGroup(parser, "Optional Options", "Those options are not mandatory for the code to work.")
    
    main.add_option("-r", "--rayfile",    help="Ray file with station pair coordinates.",                                 dest="rayfile")
    main.add_option("-o", "--output",     help="Ps Plot File",                                                            dest="psfile")
    opti.add_option("-R", "--region",     help="Region coordinates (xmin/xmax/ymin/ymax)[Default is to auto-compute].",   dest="region")
    opti.add_option("-I", "--resolution", help="Region resolution (x,y) of grids (Default is 1.0)",                       dest="resolution")
    opti.add_option("-n", "--no-clean",   help="Don't clean before recomputing.", action="store_false", default=True,     dest="clean")
    opti.add_option("-s", "--sampling",   help="Point sampling used during ray tracing over area.", default=1.0,          dest="sampling")
    opti.add_option("-m", "--mode",       help="Computation mode h: hitcount, x: cross hitcount (Default).", default='x', dest="mode")
    opti.add_option("-l", "--length",     help="Compute rays using real length.", action="store_true", default=False,     dest="length")
    opti.add_option("-u", "--uppercolor", help="Maximum of the color scale [Default is to auto-compute].", default=None,  dest="max_color")
    opti.add_option("-P", "--preview",    help="Preview the file created (Default off)", action="store_true", default=False, dest="preview")
    
    parser.add_option_group(main)
    parser.add_option_group(opti)
    
    # Parse Command Line
    #
    (options, args) = parser.parse_args()

    # Check Mode
    #
    if options.mode not in ['h', 'x']:
        parser.error("Invalid mode (%s) selected. Valid modes include one of [xh]." % options.mode)

    if options.mode == "x" and options.length: 
        parser.error("Don't know (yet) how to consider ray length in cell in cross mode.")

    # Check that there are no missing args
    #
    if args:
        parser.error("Invalid arguments: %s" % args)

    # Obtain ray file to process
    #
    if not options.rayfile:
        parser.error("No ray filename given to process.")

    if not os.path.isfile(options.rayfile):
        parser.error("Could not find the ray file indicated.")

    rayfile = os.path.basename(options.rayfile)
    raydir  = os.path.dirname(options.rayfile)

    # Read Rays
    #
    rays = parserays(options.rayfile)

    # Check ps file
    #
    if not options.psfile:
        parser.error("No output file supplied.")

    if options.psfile[-4:] == ".pdf":
        options.psfile = options.psfile[:-4] + ".ps"

    # Obtain X and Y resolution
    #
    xres = 1.0
    yres = 1.0
    if options.resolution:
        items = options.resolution.split("/")
        items = list(map(float, items))
        xres = items[0]
        yres = items[1] if len(items) > 1 else items[0]
    resolution = (xres, yres)

    # Obtain region
    #
    xmin = xmax = rays[0][0]
    ymin = ymax = rays[0][1]
    for ray in rays:
        xmin = min(xmin, ray[0], ray[2])
        xmax = max(xmax, ray[0], ray[2])
        ymin = min(ymin, ray[1], ray[3])
        ymax = max(ymax, ray[1], ray[3])

    xmin = fitnumber(xmin, xres, False)
    xmax = fitnumber(xmax, xres, True)
    ymin = fitnumber(ymin, yres, False)
    ymax = fitnumber(ymax, yres, True)

    if options.region:
        items = options.region.split("/")
        items = list(map(float, items))
        if len(items) != 4:
            parser.error("Region is invalid, needed 4 numbers.")
        xmin = items[0]
        xmax = items[1]
        ymin = items[2]
        ymax = items[3]
    limits = (xmin, xmax, ymin, ymax)

    print("Working in range -R%.3f/%.3f/%.3f/%.3f with resolution -I%f/%f" %
          (limits + resolution))

    # Building the grids from ray file
    #
    options.sampling = float(options.sampling)

    if 2 * options.sampling > min(resolution) * 111.2:
        print("*** Warning: Sampling is smaller than resolution -- aliasing will happen !!!")

    if options.length and not options.clean:
        print("*** Warning: Considering LENGTH and not recomputing grids.")

    buildgrids(rays, limits, options.sampling, resolution,
               rayfile, raydir, options.clean,
               options.length)

    pattern = os.path.join(raydir, "data-?????-%s.grd" % (rayfile))
    output  = os.path.join(raydir, "result-%s.grd" % (rayfile))
    grids   = []

    # Read grids
    #
    files = glob.glob(pattern)
    files = sorted(files)
    for f in files:
        grids.append(read(f))
        print('\rReading file: %s = %.0f%%' % (f, 100 * len(grids) / len(files)), end="")
    print()

    if len(grids) == 0:
        print('No files found for pattern "%s"' % pattern)
        sys.exit(1)

    # Print Memory Usage
    #
    print(using(" Data is Loaded"))

    # Combine rays & analyse
    #
    if options.mode == 'x':
        result = combine(grids)
    elif options.mode == 'h':
        result = hitcount(grids)

    # Decide on color scale
    #
    options.max_color = float(options.max_color) if options.max_color is not None else None
    if options.max_color is None:
        options.max_color = numpy.nanmax(result)

    # Prepare & Write result file
    #
    shutil.copy(files[0], output)
    write(output, result)

    # Print end Usage & time consumed
    #
    print(using(" Finished Processing"))

    # Plot Results
    #
    plotresults(output, options.rayfile, options.psfile, options.max_color, limits)
    options.psfile = options.psfile[:-3] + ".pdf" 
    if options.preview:
        for tool in [  "xdg-open", "evince", "atril", "okular", "gv", "gs", "open" ]:
            if os.system("which %s > /dev/null" % tool) == 0:
                print("Launching %s" % tool)
                os.system("%s %s 2> /dev/null" % (tool, options.psfile))
                break
        else:
            print("No previewing tool found.")

    # Go Out
    #
    print("Plot is done, file is %s" % options.psfile)
    sys.exit(0)

