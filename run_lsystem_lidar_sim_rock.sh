#!/bin/bash

python3 lsystem_3d_rock.py

mkdir output/7iterD1
mkdir output/7iterD2
mkdir output/7iterD3
mkdir output/7iterD4
mkdir output/7iterD5

python3 lidar_simulation_octree_rock.py
