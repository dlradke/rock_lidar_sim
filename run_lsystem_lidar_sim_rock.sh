#!/bin/bash

mkdir output/7iterVD1
mkdir output/7iterVD2
mkdir output/7iterVD3
mkdir output/7iterVD4
mkdir output/7iterVD5

python3 lsystem_3d_rock.py
python3 lidar_simulation_octree_rock.py

# git add .
# git commit -m "After 7iterD1-5 push from rock"
# git push origin master
