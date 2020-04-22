'''
    File name: lsystem_3d.py
    Author: Daniel Radke
    Date created: 30.10.2019
    Date last modified: 19.02.2020
    Python Version: 3.7.3
'''

# file to define an L-system from scratch in 3d
# functions
# h,l,u are the heading, left, up vectors of the 'turtle'

# import libraries
import random
import numpy as np
# from mayavi import mlab
import pickle
import multiprocessing as mp
import scipy


# input paramiters:
color = (105/255, 46/255, 26/255) # brown 1

# production rules
rules = {}

# bush example from Figure 1.25
global angleRange, lengthRange
iterations = 7
angleRange = [20,35]
lengthRange = [0.05,0.15]
startRadius = 0.475
axiom = 'A'
rules['A'] = "FF[&F['''^^{-f+f+f-|-f+f+f}]L!A]/////’[&F['''^^{-f+f+f-|-f+f+f}]L!A]/////’[&F['''^^{-f+f+f-|-f+f+f}]L!A]F///////’[F['''^^{-f+f+f-|-f+f+f}]L!A]"
rules['F'] = "S/////F"
rules['S'] = "F"
rules['L'] = "['''^^{-f+f+f-|-f+f+f}]"



# starting situation, origin looking up
startpos = np.asarray([0,0,0])
starth = np.asarray([0,0,1])
startl = np.asarray([0,1,0])
startu = np.asarray([-1,0,0])
starthlu = np.column_stack((starth,startl,startu))
steps = 13


# create empty ndarrays to keep pts, tris, and labels
treePts = np.empty((0,3))
treeTris = np.empty((0,3),dtype=int)
treeTrisLabel = np.empty((0,1))

# keep running total of wood volume
woodVol      = 0
woodVol1hr   = 0
woodVol10hr  = 0
woodVol100hr = 0
woodVol1000hr= 0


def rotate_turtle(axis,hlu,delta):
    # axis is either h, l, or u
    # hlu is the current hlu matrix
    if axis == 'h':
        rMat = np.asarray([[1,0,0],[0,np.cos(delta),-np.sin(delta)],[0,np.sin(delta),np.cos(delta)]])
    elif axis == 'l':
        rMat = np.asarray([[np.cos(delta),0,-np.sin(delta)],[0,1,0],[np.sin(delta),0,np.cos(delta)]])
    elif axis == 'u':
        rMat = np.asarray([[np.cos(delta),np.sin(delta),0],[-np.sin(delta),np.cos(delta),0],[0,0,1]])
    return (np.matmul(hlu,rMat))

def applyRule(input):
    output = ""
    for rule, result in rules.items():  # applying the rule by checking the current char against it
        if (input == rule):
            output = result  # Rule 1
            break
        else:
            output = input  # else ( no rule set ) output = the current char -> no rule was applied
    return output

def processString(oldStr):
    newstr = ""
    for character in oldStr:
        newstr = newstr + applyRule(character)  # build the new string
    return newstr

def createSystem(numIters, axiom):
    startString = axiom
    endString = ""
    for i in range(numIters):  # iterate with appling the rules
        print ("Iteration: {0}".format(i))
        endString = processString(startString)
        startString = endString
    return endString

def plot_cylinder(p0,p1,R,color):
    global treePts, treeTris, treeTrisLabel
    #axis and radius
    # p0 = np.array([1, 3, 2])
    # p1 = np.array([8, 5, 9])
    # R = 5
    #vector in direction of axis
    v = np.asarray(p1) - np.asarray(p0)
    #find magnitude of vector
    mag = np.linalg.norm(v)
    #unit vector in direction of axis
    v = v / mag
    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    #normalize n1
    n1 /= np.linalg.norm(n1)
    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)
    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 2)
    theta = np.linspace(0, 2 * np.pi, steps)
    # R = np.linspace(R,R*0.707,100)
    #use meshgrid to make 2d arrays
    t, theta = np.meshgrid(t, theta)
    #generate coordinates for surface
    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    # ax.plot_surface(X, Y, Z,color=color)
    # mlab.mesh(X, Y, Z, color=color)

    ## convert mesh into triangles and add to total tree triangles and points
    currNumPoints = len(treePts)
    cyPts = np.column_stack([X.flatten('F'),Y.flatten('F'),Z.flatten('F')])
    treePts = np.row_stack((treePts,cyPts))
    for p in range(0,steps-1):
        p = p + currNumPoints
        tri1 = [p,p+1,p+steps+1]
        tri2 = [p,p+steps,p+steps+1]
        treeTris = np.row_stack((treeTris,tri1))
        treeTris = np.row_stack((treeTris,tri2))
        treeTrisLabel = np.append(treeTrisLabel,('W','W'))
        # treeTrisL.append([cyPts[tri1],'W'])
        # treeTrisL.append([cyPts[tri2],'W'])

def draw_leaf(leaf):
    global treePts, treeTris, treeTrisLabel
    plotx =  np.row_stack([leaf[0:4,0],[leaf[3:7,0]]])
    ploty =  np.row_stack([leaf[0:4,1],[leaf[3:7,1]]])
    plotz =  np.row_stack([leaf[0:4,2],[leaf[3:7,2]]])
    # mlab.mesh(plotx,ploty,plotz,color=(0,1,0))

    ## convert leaf into triangles
    currNumPoints = len(treePts)
    treePts = np.row_stack((treePts,leaf))
    a = [0,1,2,3,4,5,6]
    pairs = get_consecutive_pairs1d(a)
    midLeaf = np.mean((leaf[0],leaf[3]),axis=0)
    treePts = np.row_stack((treePts,midLeaf))
    for p in pairs:
        tri = np.asarray([currNumPoints+p[0],currNumPoints+p[1],currNumPoints+7])
        treeTris = np.row_stack((treeTris,tri))
        treeTrisLabel = np.append(treeTrisLabel,('L'))


def get_rand_length():
    return np.random.uniform(lengthRange[0],lengthRange[1],1)[0]

def get_rand_angle():
    return np.radians(np.random.uniform(angleRange[0],angleRange[1],1)[0])

def get_cylinder_volume(r1,r2,h):
    r = np.mean((r1,r2))
    return np.pi * r*r * h

def drawTree(inputstring, oldpos, hlu, radius):
    # from baum.py
    global woodVol,woodVol1hr,woodVol10hr,woodVol100hr,woodVol1000hr
    i = 0 # counter for processcalculation
    processOld = 0  # old process
    newpos = oldpos
    num = []  # stack for the brackets
    # length = 1
    lastDrawR = radius
    # angleRange = [20,45]
    # angle = np.radians(22.5)

    # iterate over all characters of the string
    for character in inputstring:

        # print progress
        i += 1  # print process in percent
        process = np.floor(i * 100 / len(inputstring))
        if not process == processOld:
            print (process, "%")
            processOld = process
    # for character in string:
        # meat of the L-system. for each character, do _____
        if character == 'F':
            length = get_rand_length()
            h = hlu[:,0]
            newpos = oldpos + (length * h)
            radius = 0.97 * radius
            # mlab.plot3d([oldpos[0],newpos[0]],[oldpos[1],newpos[1]],[oldpos[2],newpos[2]])
            # ax.plot([oldpos[0],newpos[0]],[oldpos[1],newpos[1]],[oldpos[2],newpos[2]],'-',c='k',linewidth=2)
            # plt.pause(0.0000000000000000000000000001)
            if lastDrawR > radius:
                plotR = np.linspace(lastDrawR,radius,2)
            else:
                plotR = np.linspace(radius,radius,2)
            plot_cylinder(oldpos,newpos,plotR,color)
            cylVol = get_cylinder_volume(lastDrawR,radius,length)
            woodVol = woodVol + cylVol
            if lastDrawR*2 < 0.006:
                woodVol1hr = woodVol1hr + cylVol
            elif lastDrawR*2 < 0.025:
                woodVol10hr = woodVol10hr + cylVol
            elif lastDrawR*2 < 0.08:
                woodVol100hr = woodVol100hr + cylVol
            else:
                woodVol1000hr = woodVol1000hr + cylVol
            lastDrawR = radius
            oldpos = newpos

        elif character == 'f':
            h = hlu[:,0]
            newpos = oldpos + (length * h)
            oldpos = newpos
            currLeaf = np.row_stack([currLeaf,newpos])
        elif character == '+':
            # angle = np.radians(random.randrange(angleRange[0],angleRange[1]))
            # angle = a2 # page 56
            hlu = rotate_turtle('u',hlu,angle)
        elif character == '-':
            # angle = np.radians(random.randrange(angleRange[0],angleRange[1]))
            # angle = a2 # page 56
            hlu = rotate_turtle('u',hlu,-angle)
        elif character == '&':
            angle = get_rand_angle()
            # angle = a0 #page 56
            hlu = rotate_turtle('l',hlu,angle)
        elif character == '^':
            angle = get_rand_angle()
            hlu = rotate_turtle('l',hlu,-angle)
        elif character == '\\':
            angle = get_rand_angle()
            hlu = rotate_turtle('h',hlu,angle)
        elif character == '/':
            angle = get_rand_angle()
            # angle = d # page 56
            hlu = rotate_turtle('h',hlu,-angle)
        elif character == '|':
            hlu = rotate_turtle('u',hlu,np.pi)
        elif character == '!':
            # radius = 0.707 * radius
            radius = 0.75 * radius
            # radius = radius
        elif character == '[':
            num.append((oldpos, hlu, radius))
        elif character == ']':
            oldpos, hlu, radius = num.pop()
        elif character == '{':
            angle = np.radians(22.5)
            currLeaf = oldpos
            length = length * 0.8
        elif character == '}':
            draw_leaf(currLeaf)
            angle = get_rand_angle()
            length = get_rand_length()

        # # page 56
        # elif character == 'A':
        #     length = length * r1
        #     radius = radius * Wr
        # elif character == 'B':
        #     length = length * r2
        #     radius = radius * Wr
        # elif character == 'b':
        #     length = length * r1
        #     radius = radius * Wr
        # elif character == 'C':
        #     length = length * r2
        #     radius = radius * Wr
        # elif character == 'c':
        #     length = length * r1
        #     radius = radius * Wr

def get_consecutive_pairs1d(arr):
    ediff = np.ediff1d(arr)
    y = arr[:-1] + ediff
    # y = np.append(y,arr[0])
    return np.column_stack([arr[:-1],y])

def dump_pickle(data,path):
    with open(path,'wb') as filehandle:
        pickle.dump(data, filehandle)

def buildTree(treeName):
# if __name__ == '__main__':
    scipy.random.seed()
    # drawTree(createSystem(iterations, axiom), startpos)
    tree = (createSystem(iterations, axiom))
    print (tree)
    drawTree(tree, startpos, starthlu, startRadius)
    volumes = np.asarray((woodVol,woodVol1hr,woodVol10hr,woodVol100hr,woodVol1000hr))
    exportTree = [treePts,treeTris,treeTrisLabel,volumes]
    print ("Finished")
    print ("Total Wood Volume: ",woodVol)
    print ("1 hr Volume: ",woodVol1hr)
    print ("10 hr Volume: ",woodVol10hr)
    print ("100 hr Volume: ",woodVol100hr)
    print ("1000 hr Volume: ",woodVol1000hr)
    print (volumes)
    # print(treeTris.shape)
    # print(treeTris[0][0][0])
    # print(treeTris[0][0][1])
    # print(treeTris[0][0][2])
    # print(treeTris[0][1])
    ###### pickle triangles to load into ray tracer
    dump_pickle(exportTree,"data/tree_tri_export_"+str(iterations)+"iter"+str(treeName)+".data")

    # plt.show()
    # mlab.show()
    # while (1):
        # pass
        # exit()


if __name__ == '__main__':
    treeNames = ["D1","D2","D3"]
    cores = mp.cpu_count()
    pool = mp.Pool(processes=cores)
    pool.map(buildTree,treeNames)
