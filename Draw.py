import numpy as np
import pandas as pd
from PIL import Image
from mayavi import mlab
from scipy.linalg import norm
from tvtk.api import tvtk
import matplotlib.pyplot as plt

def Draw3DTreeModel(branchFileName,textureFileName):

    ## Load data
    branchDataMatrix = pd.read_csv(branchFileName,sep = ' ',header = None,names = ['ID','Point_X','Point_Y','Point_Z','CenOrder','CenRadius'])

    ## Texture
    barkImage = tvtk.JPEGReader()
    barkImage.file_name = textureFileName
    barkTexture = tvtk.Texture(input_connection = barkImage.output_port,interpolate = 1)
    resolutionH = 100
    resolutionV = 10

    ## Draw
    figSizeInch = 8
    screenDPI = 72
    mlab.figure(bgcolor = (1,1,1),size = (figSizeInch*screenDPI,figSizeInch*screenDPI))
    mlab.view(azimuth = 300,elevation = 90)
    branchIDs = np.unique(np.array(branchDataMatrix['ID']))
    for branchID in branchIDs:
        tempBranchDataMatrix = branchDataMatrix.loc[branchDataMatrix['ID'] == branchID]
        tempBranchDataMatrix_X = np.array(tempBranchDataMatrix['Point_X'])
        tempBranchDataMatrix_Y = np.array(tempBranchDataMatrix['Point_Y'])
        tempBranchDataMatrix_Z = np.array(tempBranchDataMatrix['Point_Z'])
        tempRadius = np.array(tempBranchDataMatrix['CenRadius'])
        tempRadius = tempRadius[0]
        # Initialize
        lastIndex = -1
        currentIndex = 0
        nextIndex = 1
        while nextIndex <= len(tempBranchDataMatrix_X)-1:
            # At the beginning
            if lastIndex == -1:
                currentPoint = np.array([tempBranchDataMatrix_X[currentIndex],tempBranchDataMatrix_Y[currentIndex],tempBranchDataMatrix_Z[currentIndex]])
                nextPoint = np.array([tempBranchDataMatrix_X[nextIndex],tempBranchDataMatrix_Y[nextIndex],tempBranchDataMatrix_Z[nextIndex]])
                upDirection = nextPoint - currentPoint
                for i in range(0,3):
                    if abs(upDirection[i]) < 0.000001:
                        upDirection[i] = 0.000001
                upDistance = norm(upDirection)
                upDirection = upDirection/upDistance
                # Identify X direction and Y direction
                notUpDirection = np.array([1,0,0])
                if (upDirection == notUpDirection).all():
                    notUpDirection = np.array([0,1,0])
                upXDirection = np.cross(upDirection,notUpDirection)
                upXDirection = upXDirection/norm(upXDirection)
                upYDirection = np.cross(upDirection,upXDirection)
                upYDirection = upYDirection/norm(upYDirection)
                theta = np.linspace(0,2*np.pi,resolutionH)
                bottomMesh_X = currentPoint[0] + tempRadius*np.sin(theta)*upXDirection[0] + tempRadius*np.cos(theta)*upYDirection[0]
                bottomMesh_Y = currentPoint[1] + tempRadius*np.sin(theta)*upXDirection[1] + tempRadius*np.cos(theta)*upYDirection[1]
                bottomMesh_Z = currentPoint[2] + tempRadius*np.sin(theta)*upXDirection[2] + tempRadius*np.cos(theta)*upYDirection[2]
                upMesh_X = nextPoint[0] + tempRadius*np.sin(theta)*upXDirection[0] + tempRadius*np.cos(theta)*upYDirection[0]
                upMesh_Y = nextPoint[1] + tempRadius*np.sin(theta)*upXDirection[1] + tempRadius*np.cos(theta)*upYDirection[1]
                upMesh_Z = nextPoint[2] + tempRadius*np.sin(theta)*upXDirection[2] + tempRadius*np.cos(theta)*upYDirection[2]
                for i in range(0,resolutionV):
                    if i == 0:
                        branchMesh_X = (bottomMesh_X + i*(upMesh_X-bottomMesh_X)/(resolutionV-1)).reshape(1,resolutionH)
                        branchMesh_Y = (bottomMesh_Y + i*(upMesh_Y-bottomMesh_Y)/(resolutionV-1)).reshape(1,resolutionH)
                        branchMesh_Z = (bottomMesh_Z + i*(upMesh_Z-bottomMesh_Z)/(resolutionV-1)).reshape(1,resolutionH)
                    else:
                        branchMesh_X = np.append(branchMesh_X,(bottomMesh_X + i*(upMesh_X-bottomMesh_X)/(resolutionV-1)).reshape(1,resolutionH),axis = 0)
                        branchMesh_Y = np.append(branchMesh_Y,(bottomMesh_Y + i*(upMesh_Y-bottomMesh_Y)/(resolutionV-1)).reshape(1,resolutionH),axis = 0)
                        branchMesh_Z = np.append(branchMesh_Z,(bottomMesh_Z + i*(upMesh_Z-bottomMesh_Z)/(resolutionV-1)).reshape(1,resolutionH),axis = 0)
                lastIndex = currentIndex
                currentIndex = nextIndex
                nextIndex = currentIndex + 1
            else:
                currentPoint = np.array([tempBranchDataMatrix_X[currentIndex],tempBranchDataMatrix_Y[currentIndex],tempBranchDataMatrix_Z[currentIndex]])
                bottomDirection = upDirection
                bottomXDirection = upXDirection
                bottomYDirection = upYDirection
                isContinue = False # Is there any valid next point
                while nextIndex <= len(tempBranchDataMatrix_X)-1:
                    nextPoint = np.array([tempBranchDataMatrix_X[nextIndex],tempBranchDataMatrix_Y[nextIndex],tempBranchDataMatrix_Z[nextIndex]])
                    upDirection = nextPoint - currentPoint
                    for i in range(0,3):
                        if abs(upDirection[i]) < 0.000001:
                            upDirection[i] = 0.000001
                    upDirection = upDirection/norm(upDirection)
                    crossLineDirection = np.cross(upDirection,bottomDirection)
                    crossLineDirection = crossLineDirection/norm(crossLineDirection)
                    linePoint_X = (bottomDirection[1]*np.dot(upDirection,nextPoint) - upDirection[1]*np.dot(bottomDirection,currentPoint))/(bottomDirection[1]*upDirection[0] - upDirection[1]*bottomDirection[0])
                    linePoint_Y = (np.dot(upDirection,nextPoint) - upDirection[0]*linePoint_X)/upDirection[1]
                    linePoint_Z = 0
                    linePoint = np.array([linePoint_X,linePoint_Y,linePoint_Z])
                    slopeLength = norm(currentPoint - linePoint)
                    projDistance = np.dot((currentPoint - linePoint),crossLineDirection)
                    disToLine = np.sqrt(slopeLength*slopeLength - projDistance*projDistance)
                    if disToLine < tempRadius:
                        nextIndex = nextIndex + 1
                    else:
                        isContinue = True
                        break
                if not isContinue:
                    break
                # Identify current coordinate
                notUpDirection = np.array([1,0,0])
                if (upDirection == notUpDirection).all():
                    notUpDirection = np.array([0,1,0])
                upXDirection = np.cross(upDirection,notUpDirection)
                upXDirection = upXDirection/norm(upXDirection)
                if np.dot(bottomXDirection,upXDirection) < 0:
                    upXDirection = -upXDirection
                upYDirection = np.cross(upDirection,upXDirection)
                upYDirection = upYDirection/norm(upYDirection)
                if np.dot(bottomYDirection,upYDirection) < 0:
                    upYDirection = -upYDirection
                # Calculate each slice
                theta = np.linspace(0,2*np.pi,resolutionH)
                bottomMesh_X = currentPoint[0] + tempRadius*np.sin(theta)*bottomXDirection[0] + tempRadius*np.cos(theta)*bottomYDirection[0]
                bottomMesh_Y = currentPoint[1] + tempRadius*np.sin(theta)*bottomXDirection[1] + tempRadius*np.cos(theta)*bottomYDirection[1]
                bottomMesh_Z = currentPoint[2] + tempRadius*np.sin(theta)*bottomXDirection[2] + tempRadius*np.cos(theta)*bottomYDirection[2]
                upMesh_X = nextPoint[0] + tempRadius*np.sin(theta)*upXDirection[0] + tempRadius*np.cos(theta)*upYDirection[0]
                upMesh_Y = nextPoint[1] + tempRadius*np.sin(theta)*upXDirection[1] + tempRadius*np.cos(theta)*upYDirection[1]
                upMesh_Z = nextPoint[2] + tempRadius*np.sin(theta)*upXDirection[2] + tempRadius*np.cos(theta)*upYDirection[2]
                for i in range(0,resolutionV):
                    branchMesh_X = np.append(branchMesh_X,(bottomMesh_X + i*(upMesh_X-bottomMesh_X)/(resolutionV-1)).reshape(1,resolutionH),axis = 0)
                    branchMesh_Y = np.append(branchMesh_Y,(bottomMesh_Y + i*(upMesh_Y-bottomMesh_Y)/(resolutionV-1)).reshape(1,resolutionH),axis = 0)
                    branchMesh_Z = np.append(branchMesh_Z,(bottomMesh_Z + i*(upMesh_Z-bottomMesh_Z)/(resolutionV-1)).reshape(1,resolutionH),axis = 0)
                # Update iteration condition
                lastIndex = currentIndex
                currentIndex = nextIndex
                nextIndex = currentIndex + 1
        tempMesh = mlab.mesh(branchMesh_X,branchMesh_Y,branchMesh_Z,color = (1,1,1))
        tempMesh.actor.tcoord_generator_mode = 'cylinder'
        tempMesh.actor.actor.texture = barkTexture

    ## Output
    targetDPI = 600
    outputFileName = branchFileName[:-8] + '.png'
    mlab.savefig(outputFileName,magnification = targetDPI/screenDPI) # int only
    # Change dpi without changing pixel dimensions
    outputImage = Image.open(outputFileName)
    outputImage.save(outputFileName,dpi = (targetDPI,targetDPI))
    # Interactive 3d plot
    mlab.show()

    return


if __name__ == '__main__':
    Draw3DTreeModel('./Results/Tree_1_TreeModelDraw.txt','./Data/bark.jpg')