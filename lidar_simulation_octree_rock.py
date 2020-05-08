'''
    File name: lidar_simulation_octree.py
    Author: Daniel Radke
    Date created: 25.11.2019
    Date last modified: 10.04.2020
    Python Version: 3.7.3
'''
# import libraries
# from beta_functions2 import normalize
import itertools as it
import multiprocessing as mp
import numpy as np
import os
import pickle
from pyoctree import pyoctree as ot
import time


# define functions we will need
def get_rays(openAngleD,openAngleU,openAngleL,openAngleR,nV,nH,dir):
    # heading is the direction the lidar is facing. for now it is hard coded to [1,0,0]
    # openAngle is the opening angle of the lidar scanner, in degrees
    # n is the number of rays per row/column
    windD = np.tan(np.radians(openAngleD))
    windU = np.tan(np.radians(openAngleU))
    windL = np.tan(np.radians(openAngleL))
    windR = np.tan(np.radians(openAngleR))
    windRangeV = np.linspace(-windD,windU,nV)
    windRangeH = np.linspace(-windL,windR,nH)
    if dir == 'X':
        yz = np.fromiter(it.chain(*it.product(windRangeH,windRangeV)),dtype=float).reshape(-1,2)
        rays = np.column_stack((np.repeat(1,nV*nH),yz))
    elif dir == '-X':
        yz = np.fromiter(it.chain(*it.product(windRangeH,windRangeV)),dtype=float).reshape(-1,2)
        rays = np.column_stack((np.repeat(-1,nV*nH),yz))
    elif dir == 'Y':
        xz = np.fromiter(it.chain(*it.product(windRangeH,windRangeV)),dtype=float).reshape(-1,2)
        rays = np.column_stack((xz[:,0],np.repeat(1,nV*nH),xz[:,1]))
    elif dir == '-Y':
        xz = np.fromiter(it.chain(*it.product(windRangeH,windRangeV)),dtype=float).reshape(-1,2)
        rays = np.column_stack((xz[:,0],np.repeat(-1,nV*nH),xz[:,1]))
    elif dir == 'Z':
        xy = np.fromiter(it.chain(*it.product(windRangeV,windRangeH)),dtype=float).reshape(-1,2)
        rays = np.column_stack((xy,np.repeat(1,nV*nH)))
    elif dir == '-Z':
        xy = np.fromiter(it.chain(*it.product(windRangeV,windRangeH)),dtype=float).reshape(-1,2)
        rays = np.column_stack((xy,np.repeat(-1,nV*nH)))
    rays = rays / np.linalg.norm(rays,ord=2,axis=1,keepdims=True)
    # rays = np.empty((0,3))
    # for y in windRange:
    #   for z in windRange:
    #     ray = normalize(np.asarray([1,y,z]))
    #     rays = np.row_stack((rays,ray))
    return rays

def import_pickle(path):
    with open(path,'rb') as filehandle:
        return pickle.load(filehandle)

def dump_pickle(data,path):
    with open(path,'wb') as filehandle:
        pickle.dump(data, filehandle)

def get_normal_error(sigma):
    return np.random.randn() * sigma

def rotate3(rays,delta):
    delta = np.radians(delta)
    rMat = np.asarray([[np.cos(delta),np.sin(delta),0],[-np.sin(delta),np.cos(delta),0],[0,0,1]])
    return (np.matmul(rays,rMat))

def octree_ray_intersections(origRay):
    origin = origRay[0]
    ray = origRay[1]
    # print out status
    idx = np.where((rays == ray).all(axis=1))[0][0]
    perc10 = np.round(len(rays) / 20)
    if (idx+1)%perc10 == 0:
        print("Processed Ray #: " + str(idx+1) + ". " + str(np.round((idx+1) * 100/len(rays),decimals=0)) + "%")
    beforeWoodInters = np.empty((0,2)) # array to hold intersections before a wood hit.
    rayPointList = np.row_stack((origin,origin+ray))
    rayPointList = np.asarray(rayPointList,dtype=np.float32)
    for inter in tree.rayIntersection(rayPointList):
        # beforeWoodInters = np.row_stack((beforeWoodInters,np.asarray([inter.p,inter.triLabel])))
        t = inter.s
        # tErr = t + get_normal_error(0.02) # Gaussian
        tErr = t + get_normal_error(0.00133) # Gaussian
        # tErr = t + np.random.uniform(-0.04,0.04) # Uniform
        beforeWoodInters = np.row_stack((beforeWoodInters,np.asarray([(origin + tErr*ray),inter.triLabel])))
        if (treeTrisLabel[inter.triLabel] == 'W'): # if a wood point is hit, break
            break
    return beforeWoodInters


if __name__ == '__main__':
    global tree, origin, rays, treePts, treeTris, treeTrisLabel
    # directory paths
    filePath = os.path.abspath(__file__)
    fileDir = os.path.dirname(filePath)
    # fileDir = os.getcwd()
    thesisGitDir = os.path.dirname(fileDir)

    # treeName = "7iterF3"
    for treeName in ["7iterVD1","7iterVD2","7iterVD3","7iterVD4","7iterVD5"]:
        print("Importing Tree",treeName)
        # data = import_pickle(thesisGitDir + '/data/pickle/ray_tracing/tree_tri_export_3iter.data')
        # data = import_pickle(thesisGitDir + '/l-system/MUIR/data/tree_tri_export_'+treeName+'.data')
        data = import_pickle('data/tree_tri_export_'+treeName+'.data')
        treePts = data[0]
        # treePts = treePts * 0.25
        treeTris = data[1]
        treeTrisLabel = data[2]


        # place origin of ray far enough back to hit entire trees
        xymax = max(abs(treePts[:,0].min()),abs(treePts[:,0].max()),abs(treePts[:,1].min()),abs(treePts[:,1].max()))
        # print(xymax)
        zmax = treePts[:,2].max()
        # print(zmax)
        origX = xymax + ((zmax - 1.5)/np.sqrt(3))
        # print(origX)
        '''
        # testing leaf then wood return
        treePts = np.asarray([[-1,-10,-10],[-1,-5,10],[-1,0,-10],[-1,5,10],[-1,10,-10],
                                [0,-10,-10],[0,-5,10],[0,0,-10],[0,5,10],[0,10,-10],
                                [1,-10,-10],[1,-5,10],[1,0,-10],[1,5,10],[1,10,-10]],dtype=np.float64)
        treeTris = np.asarray([[0,1,2],[1,2,3],[2,3,4],
                                [5,6,7],[6,7,8],[7,8,9],
                                [10,11,12],[11,12,13],[12,13,14]])
        treeTrisLabel = np.asarray([['L'],['L'],['L'],['L'],['W'],['L'],['W'],['W'],['W']])
        treeTrisLabel = treeTrisLabel.reshape(-1,)
        '''
        t = time.time()
        print("Building Octree")
        print("Number of triangle: " + str(len(treeTris)))
        tree = ot.PyOctree(treePts,np.asarray(treeTris,dtype=np.int32))

        for i in ([[1,1001,901],[2,2001,1801]]):
            k = i[0]
            nV = i[1]
            nH = i[2]
            # nV = 2001
            # nH = 1801
            origin1 = np.asarray([-origX,0,1.5])
            origin2 = rotate3(origin1,120)
            origin3 = rotate3(origin1,240)
            print("Creating Rays")
            print("Position of First Lidar Scanner: ",origin1)
            rays1 = get_rays(40,60,45,45,nV,nH,'X')
            rays2 = rotate3(rays1,120)
            rays3 = rotate3(rays1,240)
            rays = np.row_stack((rays1,rays2,rays3))
            origins = np.row_stack((np.tile(origin1,(len(rays1),1)),np.tile(origin2,(len(rays2),1)),np.tile(origin3,(len(rays3),1))))
            origRays = np.column_stack((origins,rays)).reshape(-1,2,3)



            print("Number of rays: " + str(len(rays)))
            # parallelize
            cores = mp.cpu_count()
            pool = mp.Pool(processes=cores)
            print("Starting Ray Tracing")


            intersectionPts = pool.map(octree_ray_intersections,origRays)

            print("Finished Ray Tracing")
            elapsed = time.time() - t
            print('Elapsed Time: ' + str(elapsed) + ' seconds')
            print('Elapsed Time: ' + str(elapsed/60) + ' minutes')
            pool.terminate()

            # formate points and cooresponding triangles
            intersectionPts = np.row_stack(intersectionPts)
            intersectionTris = np.asarray(intersectionPts[:,1],dtype=np.int32)
            intersectionPts = np.row_stack(intersectionPts[:,0])
            print("Number of intersections: " + str(len(intersectionPts)))

            # seperate into wood and leaf points
            woodPts = intersectionPts[treeTrisLabel[intersectionTris] == 'W']
            leafPts = intersectionPts[treeTrisLabel[intersectionTris] == 'L']

            # dump pickle
            dumpPath = 'output/'+treeName+'/'+str(nV)+'x'+str(nH)+'rays.data'
            dump_pickle([np.asarray(intersectionPts,dtype='float32'),treeTrisLabel[intersectionTris],origin1],dumpPath)
            np.save('output/'+treeName+'/'+str(k)+'k/coords.npy',np.asarray(intersectionPts,dtype='float32'))


    # from playsound import playsound
    # playsound('/Users/dlradke/Documents/misc/sounds/sncf_2005.m4r')
    # plot
    # '''
    # from mayavi import mlab
    # mlab.points3d(treePts[:,0],treePts[:,1],treePts[:,2],mode='point')
    # mlab.points3d(woodPts[:,0],woodPts[:,1],woodPts[:,2],color=(105/255, 46/255, 26/255),mode='sphere',scale_factor=0.01)
    # mlab.points3d(leafPts[:,0],leafPts[:,1],leafPts[:,2],color=(0,1,0),mode='sphere',scale_factor=0.01)
    # mlab.points3d(woodPts[:,0],woodPts[:,1],woodPts[:,2],color=(105/255, 46/255, 26/255),mode='point')
    # mlab.points3d(leafPts[:,0],leafPts[:,1],leafPts[:,2],color=(0,1,0),mode='point')
    # mlab.show()
    '''

    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.plot(woodPts[:,0],woodPts[:,1],woodPts[:,2],'o',color='r')
    plt.plot(leafPts[:,0],leafPts[:,1],leafPts[:,2],'o',color='g')
    plt.show()
    # '''
