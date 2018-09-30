# Que Raios

Project to develop a ray cross indicator (RCI) quality indicator for
evaluating ANT maps.

The code reads a 4 column file with coordinates (lon, lat, lon, lat) of pair
of stations that were used to perform ANT.  It computes the RCI indicator
for the pairs and plots a Qc map.  Computations are done in GMT/netcdf
compatible grids.  It depends on GMT for plotting and doing some initial
grid computation.

## Distribution

This code is GPL. It can be obtained from: https://github.com/marcelobianchi/queraios

## Usage

The code has a built in help message that can be used to get a quick start.
Just use *python queraios.py*. The input file is a 4-column file containing
the stations pair coordinates.

```
Usage: queraios.py [options]

Options:
  -h, --help            show this help message and exit

  Main Mandatory Options:
    Those options are mandatory for the code to work.

    -r RAYFILE, --rayfile=RAYFILE
                        Ray file with station pair coordinates.
    -o PSFILE, --output=PSFILE
                        Ps Plot File

  Optional Options:
    Those options are not mandatory for the code to work.

    -R REGION, --region=REGION
                        Region coordinates (xmin/xmax/ymin/ymax)[Default is to
                        auto-compute].
    -I RESOLUTION, --resolution=RESOLUTION
                        X/Y cell resolution in degrees (Default is -I1.0/1.0)
    -n, --no-clean      Don't clean before recomputing.
    -s SAMPLING, --sampling=SAMPLING
                        Point sampling used during ray tracing over area (km).
    -m MODE, --mode=MODE
                        Computation mode h: hitcount, x: cross hitcount
                        (Default).
    -l, --length        Compute rays using real length. Normally a ray passing
                        a cell counts 1, if this option is given, the ray will
                        count the amount of km it was inside the cell.
    -u MAX_COLOR, --uppercolor=MAX_COLOR
                        Maximum of the color scale [Default is to auto-
                        compute].
    -P, --preview       Preview the file created (Default off)
```
