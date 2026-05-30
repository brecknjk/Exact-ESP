# Continuous Dijkstra with Pruning for ESP
This project can extracts coastlines from osm-data in a pbf-format and create constrained Delaunay Triangulation, that can be used to route ships around the globe.
A form of a continuous Dijkstra was implemented inspired by Polyanya (https://www.ijcai.org/proceedings/2017/0070.pdf).

## Content of this repo
- settings.json
Because this project has many configurable parts, a settings file is used for controlling all relevant parameters.
Some parameters can be overwritten via command line arguments. For example the verbosity level, the debug level and importantly the path to the settings.json file.
- scrips/
In this folder are some small helper scrips used to visualize and a shrinker that reduces the complexity of the obstacles.
- data/
Here is the input data saved and the default location for the produced output data.
To keep the size of this project reasonable small only smaller instances are saved in git.
If you need access to all instances please contact us.
- cpp_program/extractor
This program creates and outputs a constrained delaunay triangulation created from the obstacle polylines.
If USE_OSMPBF is set its possible to parse a pbf-file in order to extract the coastlines for OSM-data.
With the CLI arg "-wcPa" it is possible to create a worst-case scenario for polyanya for arbitrary sizes.
- cpp_program/router
This is the main program of this project.

## Setup
```
git clone --recursive https://github.tik.uni-stuttgart.de/st167480/ContinuousDijkstra.git
mkdir build
cd build
cmake -DCGAL_USE_EPECK=false -DCGAL_DIR="~/src/cgal" -DCMAKE_BUILD_TYPE=Release ../
make
```
The "CGAL_USE_EPECK" flag is used to decide if exact constrution and sqrt should be used. 
If set to true this will greatly increase the runtime, but the calculated routs are exact.
If cmake couldn't find the CGAL folder it can be specified via CGAL_DIR.

### how to use
```
./build/src/router/router ./data/triangulation/coastlines_shrunk_2.fmi benchmark -m -tr 100 -d 0 -v 1 -fp ./
```

```
./build/src/router/router ./data/triangulation/aurora.mesh batch -tr 100 -d 0 -v 1 -fp ./
```

## Router
### Input File Style
```
nodeCount
EdgeCount
n.lat n.lon
...
e.n0 e.n1 e.e0' e.e1' // it holds that e[e.e0' / 2].n[e.e0' mod 2] == e.n0

```

## Speedups
- A*-Heuristik funktion to guide the search in two versions a better but harder to compute and a simpler one.
- Deadend pruning: e=(a,b) is deadend iff a and b correspond to the same obstical -> don't look inside if goal isn't
- duplicate detection (aka if cost[generator] has changed)
- simple redudant pruning: for e=(a,b) if cost[a] + dist(a,b) < cost[generator] + dist(generator,b) -> prune and vice versa
- just generate from a point with the same cost once

## Running the original polyanya implementation
'''
git clone https://bitbucket.org/dharabor/pathfinding.git
cd ./anyangle/polyanya
make
./bin/scenariorunner ./meshes/arena.mesh ./scenarios/arena.scen
'''
scen:
version 1
0	maps/dao/arena.map	49	49	1	11	1	12	1
0	maps/dao/arena.map	49	49	1	12	1	10	2 // s.bucket >> map >> s.xsize >> s.ysize >> s.start.x >> s.start.y >> s.goal.x >> s.goal.y >> s.gridcost
... 

mesh:
"""
mesh
2
numberNodes numberPolygons
x y neighbourCount p0 p1 ... // points
cornerCount n0 n1 ... p0 p1 ...                   //polygons n0 n1 -> is edge to p0

2 2 2 -1 4
1 3 2 15 -1
2 3 6 1 15 -1 4 0 14
3 2 5 2 0 4 -1 5
3 1 2 5 -1
15 1 3 -1 2 5
15 3 6 6 3 14 0 2 -1
...
3 60 23 91 96 115 93
3 61 62 63 119 116 -1
3 61 63 76 82 118 84
"""