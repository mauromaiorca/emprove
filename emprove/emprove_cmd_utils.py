#!/usr/bin/python3


import argparse
import os.path
import numpy as np
import emprove_core
from emprove import starHandler
from os import PathLike,makedirs
from emprove import utils
from numpy.fft import fftn, ifftn
import pandas as pd

emprove_parser = argparse.ArgumentParser(
    prog="emprove_utils",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
command = emprove_parser.add_subparsers(dest="command")


#################################
## emprove_maskedCrop 
emprove_maskedCrop = command.add_parser (
    "maskedCrop", description="compute automatic crop of an image based on mask", help='compute automatic crop of an image based on mask'
)
emprove_maskedCrop.add_argument("--i", required=True, type=str, help="file with the input mrc map")
emprove_maskedCrop.add_argument("--mask", required=True, type=str, default="", help="file with the input mask mrc map")
emprove_maskedCrop.add_argument("--padding", required=False, type=int, default=2, help="padding value for the mask, default=2")
emprove_maskedCrop.add_argument("--o", required=True, type=str, help="outputFilename")
def maskedCrop(args):
    sizeMap=emprove_core.sizeMRC(args.i)
    sizeMask=emprove_core.sizeMRC(args.mask)
    inputMap=np.array(emprove_core.ReadMRC(args.i))
    inputMap=np.reshape(inputMap,sizeMap)
    inputMask=np.array(emprove_core.ReadMRC(args.mask))
    inputMask=np.reshape(inputMask,sizeMap)

    nonzero_indices = np.where(inputMask > 0.1)
    min_z = np.min(nonzero_indices[0])
    max_z = np.max(nonzero_indices[0]) + 1
    min_y = np.min(nonzero_indices[1])
    max_y = np.max(nonzero_indices[1]) + 1
    min_x = np.min(nonzero_indices[2])
    max_x = np.max(nonzero_indices[2]) + 1

    max_size = max(max_x - min_x, max_y - min_y, max_z - min_z)

    # Adjust bounding box
    half_size = max_size // 2
    center_z = (max_z + min_z) // 2
    center_y = (max_y + min_y) // 2
    center_x = (max_x + min_x) // 2
    min_z = center_z - half_size
    max_z = center_z + half_size + max_size % 2
    min_y = center_y - half_size
    max_y = center_y + half_size + max_size % 2
    min_x = center_x - half_size
    max_x = center_x + half_size + max_size % 2

    # Introduce padding
    padding = args.padding
    min_z -= padding
    max_z += padding
    min_y -= padding
    max_y += padding
    min_x -= padding
    max_x += padding

    # Ensure the bounding box doesn't go outside the original image
    min_z = max(min_z, 0)
    max_z = min(max_z, inputMap.shape[0])
    min_y = max(min_y, 0)
    max_y = min(max_y, inputMap.shape[1])
    min_x = max(min_x, 0)
    max_x = min(max_x, inputMap.shape[2])
    croppedMap = inputMap[min_z:max_z, min_y:max_y, min_x:max_x]

    print("Size:", inputMap.size)
    print("Shape:", inputMap.shape)
    print("Data Type:", inputMap.dtype)
    emprove_core.WriteMRC(croppedMap.flatten().tolist(), args.o ,croppedMap.shape[2],croppedMap.shape[1],croppedMap.shape[0],1)

    #print("sizeMap=",sizeMap)
    #emprove_core.WriteEmptyMRC(outputStackBasename+'.mrcs',sizeI[0],sizeI[1],len(imageNames[imageNameTag]))



#################################
## rotation_average_2d 
emprove_rotation_average_2d = command.add_parser (
    "rotation_average_2d", description="rotation_average_2d", help='rotation_average_2d'
)
emprove_rotation_average_2d.add_argument("--i", required=True, type=str, help="2D input map")
emprove_rotation_average_2d.add_argument("--csv", required=False, type=str, default="", help="output csv file with results")
emprove_rotation_average_2d.add_argument("--o", required=False, type=str, default="",  help="output 1D mrc map with radial average")
def rotation_average_2d(args):
    import csv
    sizeMap=emprove_core.sizeMRC(args.i)
    inputMap=np.array(emprove_core.ReadMRC(args.i))
    inputMap=np.reshape(inputMap,sizeMap)
    center = np.array([inputMap.shape[0]//2, inputMap.shape[1]//2])
    y, x = np.ogrid[:inputMap.shape[0], :inputMap.shape[1]]
    distances_squared = (x - center[1])**2 + (y - center[0])**2
    distances = np.sqrt(distances_squared).astype(np.int64)
    radial_sum = np.bincount(distances.ravel(), weights=inputMap.ravel())
    radial_count = np.bincount(distances.ravel())
    radial_count[radial_count == 0] = 1
    radial_avg = radial_sum / radial_count
    # Plot the radial average
    if not args.o == "":
        emprove_core.WriteMRC(radial_avg.flatten().tolist(), args.o ,len(radial_avg),1,1,1)
    elif not args.csv == "":
        with open(args.csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)  # Pass the opened file object, not a string
            writer.writerow(radial_avg)
    else:
        import matplotlib.pyplot as plt
        plt.plot(radial_avg)
        plt.xlabel('Distance from Center')
        plt.ylabel('Average Intensity')
        plt.title('Radial Average')
        plt.grid(True)
        plt.show()





#################################
## emprove_rotation_average_stack 
emprove_rotation_average_stack = command.add_parser (
    "rotation_average_stack", description="rotation_average_stack", help='rotation_average_stack'
)
emprove_rotation_average_stack.add_argument("--i", required=True, type=str, help="2D stack")
emprove_rotation_average_stack.add_argument("--csv", required=True, type=str, default="", help="output csv file with results")
emprove_rotation_average_stack.add_argument("--o", required=False, type=str, default="",  help="output 1D mrc map with radial average")
# Plotting function
def plot_rotation_averageStack(filename, method='euclidean'):
    import csv
    import matplotlib.pyplot as plt
    curves = []
    with open(filename, 'r', newline='') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            curves.append([float(val) for val in row])
    curves = np.array(curves)
    avg_curve = np.mean(curves, axis=0)
    # Compute pairwise differences
    max_distance = float('-inf')
    most_diff_curves = (None, None)
    for i in range(len(curves)):
        for j in range(i+1, len(curves)):
            if method == 'euclidean':
                distance = np.linalg.norm(curves[i] - curves[j])
            elif method == 'dtw':
                import fastdtw
                distance, _ = fastdtw.fastdtw(curves[i], curves[j])
            elif method == 'cosine':
                distance = 1 - np.dot(curves[i], curves[j]) / (np.linalg.norm(curves[i]) * np.linalg.norm(curves[j]))

            if distance > max_distance:
                max_distance = distance
                most_diff_curves = (curves[i], curves[j])
    x_vals = list(range(len(avg_curve)))
    plt.plot(x_vals, most_diff_curves[0], label="Most Different 1")
    plt.plot(x_vals, most_diff_curves[1], label="Most Different 2")
    plt.plot(x_vals, avg_curve, label="Average Curve", linestyle='--')
    plt.xlabel('Distance from Center')
    plt.ylabel('Average Intensity')
    plt.title(f'Radial Averages using {method} method')
    plt.legend()
    plt.show()
def rotation_average_stack(args):
    import csv
    if os.path.isfile(args.csv):
        os.remove(args.csv)
    sizeMap=emprove_core.sizeMRC(args.i)
    print ("size=", sizeMap[2])
    #emprove_core.WriteEmptyMRC(outputStackBasename+'.mrcs',sizeI[0],sizeI[1],len(imageNames[imageNameTag]))
    for ii in range(0,sizeMap[2]):
        #print (ii)
        inputMap=emprove_core.ReadMrcSlice(args.i,ii)
        inputMap=np.reshape(inputMap,(sizeMap[0],sizeMap[1]))
        center = np.array([inputMap.shape[0]//2, inputMap.shape[1]//2])
        y, x = np.ogrid[:inputMap.shape[0], :inputMap.shape[1]]
        distances_squared = (x - center[1])**2 + (y - center[0])**2
        distances = np.sqrt(distances_squared).astype(np.int64)
        radial_sum = np.bincount(distances.ravel(), weights=inputMap.ravel())
        radial_count = np.bincount(distances.ravel())
        radial_count[radial_count == 0] = 1
        radial_avg = radial_sum / radial_count
        # Plot the radial average
        with open(args.csv, 'a', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(radial_avg)
    plot_rotation_averageStack(args.csv, method='dtw')





emprove_project_volume = command.add_parser (
    "project_volume", description="project a volume, based on the star file", help='project a volume, based on the star file'
)
emprove_project_volume.add_argument("--star", required=False, type=str, help="input star file")
emprove_project_volume.add_argument("--map", required=True, type=str, help="file with the input mrc map")
emprove_project_volume.add_argument("--rot", required=False, type=float, default=0, help="rot angle")
emprove_project_volume.add_argument("--tilt", required=False, type=float, default=0, help="tilt angle")
emprove_project_volume.add_argument("--psi", required=False, type=float, default=0, help="psi angle")
emprove_project_volume.add_argument("--xoff", required=False, type=float, default=0, help="xoff in pixels")
emprove_project_volume.add_argument("--yoff", required=False, type=float, default=0, help="yoff in pixels")
emprove_project_volume.add_argument("--c_style", action="store_true", help="c star like projection")
emprove_project_volume.add_argument("--o", required=True, type=str, help="output for mrcs file")
def project_volume(args):
    print("project_volume")
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")  # use GPU 0 if available, else CPU

    sizeMap=emprove_core.sizeMRC(args.map)
    inputMap=np.array(emprove_core.ReadMRC(args.map))
    inputMap=np.reshape(inputMap,(sizeMap[2],sizeMap[1],sizeMap[0]))
    rlnAngleRot= torch.tensor(np.radians(args.rot), dtype=torch.float32).to(device)
    rlnAngleTilt= torch.tensor(np.radians(args.tilt), dtype=torch.float32).to(device)
    rlnAnglePsi= torch.tensor(np.radians(args.psi), dtype=torch.float32).to(device)
    inputMapTorch = torch.tensor(inputMap, dtype=torch.float32)

    if not args.star == "":
        print ("star file") 
        version=starHandler.infoStarFile(args.star)[2]
        suffix_particles='particles'
        if not version=="relion_v31":
            suffix_particles=''
        image_names = starHandler.read_star_columns_from_sections(args.star, suffix_particles, '_rlnImageName')
        numElements=image_names.shape[0]
        emprove_core.WriteEmptyMRC(args.o,sizeMap[0],sizeMap[1],numElements)

        if version=="relion_v31":
            print ("relion_v31")
            print ("file=",args.star)
            headerData = starHandler.read_star_columns_from_sections(args.star, 'optics', '_rlnImagePixelSize')
            #print ("headerData=",headerData)
            particle_parameters = starHandler.read_star_columns_from_sections(args.star, 'particles', ['_rlnAngleRot','_rlnAngleTilt', '_rlnAnglePsi','_rlnOriginXAngst', '_rlnOriginYAngst'])
            #print ("particle_parameters=",particle_parameters)

            for index, row in particle_parameters.iterrows():
                if args.c_style:
                    print('c-style projections')
                    phi = float(row['_rlnAngleRot'])
                    theta = float(row['_rlnAngleTilt']) 
                    psi = float(row['_rlnAnglePsi']) 
                    projection=emprove_core.projectMap(inputMap.flatten().tolist(),sizeMap[0],sizeMap[1],sizeMap[2],float(row['_rlnAngleRot']),float(row['_rlnAngleTilt']),float(row['_rlnAnglePsi']),float(row['_rlnOriginXAngst']), float(row['_rlnOriginYAngst']),0)
                    emprove_core.ReplaceMrcSlice(projection.flatten().tolist(),args.o,sizeMap[0],sizeMap[1],index)
                else:
                    phi = torch.tensor(np.radians(float(row['_rlnAngleRot'])), dtype=torch.float32).to(device)
                    theta = torch.tensor(np.radians(float(row['_rlnAngleTilt'])), dtype=torch.float32).to(device)
                    psi = torch.tensor(np.radians(float(row['_rlnAnglePsi'])), dtype=torch.float32).to(device)
                    projection = projector_torch.forward_projection(inputMapTorch, phi, theta, psi, float(row['_rlnOriginXAngst']), float(row['_rlnOriginYAngst']), device=device)
                    projectionOut = projection.cpu().numpy()
                    emprove_core.ReplaceMrcSlice(projectionOut,args.o,sizeMap[0],sizeMap[1],index)
        else:
            particle_parameters = starHandler.read_star_columns_from_sections(args.star, '', ['_rlnAngleRot','_rlnAngleTilt', '_rlnAnglePsi','_rlnOriginX', '_rlnOriginY' ])
            #print ("particle_parameters=",particle_parameters)
            for index, row in particle_parameters.iterrows():
                if args.c_style:
                    print('c-style projections')
                    phi = float(row['_rlnAngleRot'])
                    theta = float(row['_rlnAngleTilt']) 
                    psi = float(row['_rlnAnglePsi']) 
                    projection=emprove_core.projectMap(inputMap.flatten().tolist(),sizeMap[0],sizeMap[1],sizeMap[2],float(row['_rlnAngleRot']),float(row['_rlnAngleTilt']),float(row['_rlnAnglePsi']),float(row['_rlnOriginX']), float(row['_rlnOriginY']),0)
                    emprove_core.ReplaceMrcSlice(projection,args.o,sizeMap[0],sizeMap[1],index)
                else:
                    phi = torch.tensor(np.radians(float(row['_rlnAngleRot'])), dtype=torch.float32).to(device)
                    theta = torch.tensor(np.radians(float(row['_rlnAngleTilt'])), dtype=torch.float32).to(device)
                    psi = torch.tensor(np.radians(float(row['_rlnAnglePsi'])), dtype=torch.float32).to(device)
                    projection = projector_torch.forward_projection(inputMapTorch, phi, theta, psi, float(row['_rlnOriginX']), float(row['_rlnOriginY']), device=device)
                    projectionOut = projection.cpu().numpy()
                    emprove_core.ReplaceMrcSlice(projectionOut.flatten().tolist(),args.o,sizeMap[0],sizeMap[1],index)

    else:
        projection = projector_torch.forward_projection(inputMapTorch, rlnAngleRot, rlnAngleTilt, rlnAnglePsi, args.xoff, args.yoff, device=device)
        projectionOut = projection.cpu().numpy()
        emprove_core.WriteMRC(projection.flatten().tolist(), args.o ,projectionOut.shape[1],projectionOut.shape[0],1,1)




def backward_projection(data, angle, device):
    # Create an empty volume for the backprojection
    volume_size = data.shape[0] # assuming square volume for simplicity
    backproj = torch.zeros((volume_size, volume_size), device=device)
    
    # Calculate the projection for each pixel
    for x in range(volume_size):
        for y in range(volume_size):
            proj_x = int(x * torch.cos(torch.radians(angle)) - y * torch.sin(torch.radians(angle)))
            if 0 <= proj_x < volume_size:
                backproj[x, y] = data[proj_x]
                
    return backproj




def apply_2d_ramp_filter(projection_fft, nx, ny, device='cpu'):
    """
    Applies a 2D ramp filter to the provided FFT of a projection.

    Parameters:
    - projection_fft (torch.Tensor): 2D Fourier transform of the projection.
    - nx (int): The size of the projection in the x-direction.
    - ny (int): The size of the projection in the y-direction.
    - device (str, optional): Device to be used for computations.

    Returns:
    - torch.Tensor: The filtered projection in the frequency domain.
    """
    
    # 1. Create 1D ramp filters
    ramp_1d_x = torch.fft.fftfreq(nx, device=device).abs()
    ramp_1d_y = torch.fft.fftfreq(ny, device=device).abs()
    
    # 2. Create 2D ramp filter from the 1D ramp filters
    ramp_2d = torch.outer(ramp_1d_y, ramp_1d_x)
    
    # Move ramp_2d to the desired device
    ramp_2d = ramp_2d.to(device)

    #print("projection_fft shape:", projection_fft.shape)
    #print("ramp_2d shape:", ramp_2d.shape)

    # 3. Multiply the 2D FFT of the projection by the 2D ramp filter
    filtered_projection_fft = projection_fft.reshape(ny, nx) * ramp_2d

    # 4. Return the filtered FFT of the projection
    return filtered_projection_fft


def euler_to_rotation_matrix(phi, theta, psi, device='cpu'):
    """
    Returns rotation matrix from Euler angles.
    """
    R_x = torch.tensor([[1, 0, 0],
                        [0, torch.cos(phi), -torch.sin(phi)],
                        [0, torch.sin(phi), torch.cos(phi)]], device=device)
    R_y = torch.tensor([[torch.cos(theta), 0, torch.sin(theta)],
                        [0, 1, 0],
                        [-torch.sin(theta), 0, torch.cos(theta)]], device=device)
    R_z = torch.tensor([[torch.cos(psi), -torch.sin(psi), 0],
                        [torch.sin(psi), torch.cos(psi), 0],
                        [0, 0, 1]], device=device)
    R = R_z @ R_y @ R_x
    return R



def fftshift_1d_torch(tensor, dim=0):
    N = tensor.size(dim)
    split_tensors = torch.split(tensor, [N//2, N - N//2], dim)
    return torch.cat((split_tensors[1], split_tensors[0]), dim=dim)



def FBP_reconstructionFFT(starFile, outputMRC, device='cpu'):
    # Presuming some missing parts of the code based on context
    version=starHandler.infoStarFile(starFile)[2]
    print("star file version=",version)
    
    if version=="relion_v31":
        suffix_particles='particles'
        particle_parameters = starHandler.read_star_columns_from_sections(starFile, suffix_particles,  ['_rlnImageName','_rlnAngleRot','_rlnAngleTilt',  '_rlnAnglePsi','_rlnOriginXAngst', '_rlnOriginYAngst'])
    else:
        suffix_particles=''
        particle_parameters = starHandler.read_star_columns_from_sections(starFile, suffix_particles,  ['_rlnImageName','_rlnAngleRot','_rlnAngleTilt',  '_rlnAnglePsi','_rlnOriginX', '_rlnOriginY'])

    image_name_first_row = str(particle_parameters['_rlnImageName'].iloc[0])
    firstStackName = image_name_first_row[image_name_first_row.find('@')+1:]
    sizeMap = emprove_core.sizeMRC(firstStackName)
    nx, ny, nz = sizeMap[0], sizeMap[1], max(sizeMap[0], sizeMap[1])
    
    outmap=torch.zeros((nx, ny, nz), device=device)
    outmap_complex = torch.complex(outmap, torch.zeros_like(outmap, device=device))
    outmap_fft = torch.fft.fftn(torch.fft.fftshift(outmap_complex))
    countmap = torch.zeros_like(outmap_fft, device=device).real
    countmap = torch.zeros_like(outmap_fft, device=device).real

    # Create a grid for the 3D coordinates
    x = torch.linspace(0, nx-1, nx, device=device)
    y = torch.linspace(0, ny-1, ny, device=device)
    z = torch.linspace(0, nz-1, nz, device=device)
    xx, yy, zz = torch.meshgrid(x, y, z, indexing='xy')
    # Translate to have center at (0, 0, 0)
    xx = xx - nx/2
    yy = yy - ny/2
    zz = zz - nz/2

    for index, row in particle_parameters.iterrows():
        tmpLine = row['_rlnImageName']
        atPosition = tmpLine.find('@')
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition+1:]
        phi = torch.tensor(np.radians(float(row['_rlnAngleRot']) ), dtype=torch.float32).to(device)
        theta = torch.tensor(np.radians(float(row['_rlnAngleTilt']) ), dtype=torch.float32).to(device)
        psi = torch.tensor(np.radians(float(row['_rlnAnglePsi']) ), dtype=torch.float32).to(device)
        if version=="relion_v31":
            shiftX = float(row['_rlnOriginXAngst']) 
            shiftY = float(row['_rlnOriginYAngst'])
        else:
            shiftX = float(row['_rlnOriginX']) 
            shiftY = float(row['_rlnOriginY'])

        p = np.array(emprove_core.ReadMrcSlice(stackName, imageNo-1), dtype=np.float32)

        p_torch=torch.tensor(p).to(device)
        p_torch_fft=torch.fft.fftshift(torch.fft.fftn(p_torch))
        #p_torch_fft=apply_2d_ramp_filter(p_torch_fft, nx, ny, device=device)
        p_torch_fft=fftshift_1d_torch(p_torch_fft,0)
        #p_torch_ifft=torch.fft.ifftn(torch.fft.ifftshift(p_torch_fft)).real
        #emprove_core.WriteMRC(p_torch_ifft.cpu().numpy().flatten().tolist(), "ciao1.mrc" ,nx, ny, 1,1)


        # Create a grid for the 2D coordinates of p_torch
        x2d = torch.linspace(0, nx-1, nx, device=device)
        y2d = torch.linspace(0, ny-1, ny, device=device)
        xx2d, yy2d = torch.meshgrid(x2d, y2d, indexing='xy')
        # Translate to have center at (0, 0)
        xx2d = xx2d - nx/2
        yy2d = yy2d - ny/2        
        # Flatten and stack 2D coordinates, assume z=0 for 2D
        coords2d = torch.stack([xx2d.flatten(), yy2d.flatten(), torch.zeros_like(xx2d).flatten(), torch.ones_like(xx2d).flatten()])
        # Apply rotation
        R = compute_affine_matrix(phi, theta, psi, shiftX, shiftY, invertMatrix=True, device=device)

        rotated_coords2d = R @ coords2d
        # Map rotated coordinates to volume coordinates
        rotated_coords2d[0] += nx/2
        rotated_coords2d[1] += ny/2
        rotated_coords2d[2] += nz/2
        # Create mask based on the rotated 2D coordinates
        mask2d = (rotated_coords2d[0] >= 0) & (rotated_coords2d[0] < nx) & (rotated_coords2d[1] >= 0) & (rotated_coords2d[1] < ny) &  (rotated_coords2d[2] >= 0) & (rotated_coords2d[2] < nz)
        valid_indices2d = rotated_coords2d[:, mask2d]
        p_complex_values = p_torch_fft.flatten()[mask2d]
        outmap_fft[valid_indices2d[0].long(), valid_indices2d[1].long(), valid_indices2d[2].long()] += p_complex_values
        countmap[valid_indices2d[0].long(), valid_indices2d[1].long(), valid_indices2d[2].long()] += 1
        p_real_values = p_torch.flatten()[mask2d]
        outmap[valid_indices2d[0].long(), valid_indices2d[1].long(), valid_indices2d[2].long()] += p_real_values

    mask_non_zero = countmap > 0
    outmap[mask_non_zero] /= countmap[mask_non_zero]

    outmap=torch.fft.ifftn(outmap_fft).real
    emprove_core.WriteMRC(outmap.cpu().numpy().flatten().tolist(), outputMRC ,nx, ny, nz,1)


def FBP_reconstruction(starFile, outputMRC, device='cpu'):
    # Presuming some missing parts of the code based on context
    version=starHandler.infoStarFile(starFile)[2]
    print("star file version=",version)
    
    if version=="relion_v31":
        suffix_particles='particles'
        particle_parameters = starHandler.read_star_columns_from_sections(starFile, suffix_particles,  ['_rlnImageName','_rlnAngleRot','_rlnAngleTilt',  '_rlnAnglePsi','_rlnOriginXAngst', '_rlnOriginYAngst'])
    else:
        suffix_particles=''
        particle_parameters = starHandler.read_star_columns_from_sections(starFile, suffix_particles,  ['_rlnImageName','_rlnAngleRot','_rlnAngleTilt',  '_rlnAnglePsi','_rlnOriginX', '_rlnOriginY'])

    image_name_first_row = str(particle_parameters['_rlnImageName'].iloc[0])
    firstStackName = image_name_first_row[image_name_first_row.find('@')+1:]
    sizeMap = emprove_core.sizeMRC(firstStackName)
    nx, ny, nz = sizeMap[0], sizeMap[1], max(sizeMap[0], sizeMap[1])
    
    outmap=torch.zeros((nx, ny, nz), device=device)
    for index, row in particle_parameters.iterrows():
        tmpLine = row['_rlnImageName']
        atPosition = tmpLine.find('@')
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition+1:]
        phi = torch.tensor(np.radians(float(row['_rlnAngleRot']) ), dtype=torch.float32).to(device)
        theta = -torch.tensor(np.radians(float(row['_rlnAngleTilt']) ), dtype=torch.float32).to(device)
        psi = torch.tensor(np.radians(float(row['_rlnAnglePsi']) ), dtype=torch.float32).to(device)
        if version=="relion_v31":
            shiftX = float(row['_rlnOriginXAngst']) 
            shiftY = float(row['_rlnOriginYAngst'])
        else:
            shiftX = float(row['_rlnOriginX']) 
            shiftY = float(row['_rlnOriginY'])

        # Load the 2D projection as a tensor
        p=np.array(emprove_core.ReadMrcSlice(stackName, imageNo-1))
        p=np.reshape(p,(nx,ny))
        p_torch = torch.tensor(p, device=device)
        # Create a grid for the 2D coordinates of p
        x2d = torch.linspace(-nx/2, nx/2-1, nx, device=device)
        y2d = torch.linspace(-ny/2, ny/2-1, ny, device=device)
        xx2d, yy2d = torch.meshgrid(x2d, y2d, indexing='xy')
        # Translate to have center at (0, 0)
        xx2d = xx2d - nx/2
        yy2d = yy2d - ny/2
        
        # Create 3D coordinates for backprojection (essentially smearing the 2D image into 3D space)
        xx3d, yy3d, zz3d = torch.meshgrid(x2d, y2d, torch.linspace(-nz/2, nz/2-1, nz, device=device), indexing='xy')
        # Clone the tensors before translating to avoid in-place operations on view tensors
        xx3d = xx3d.clone() + nx/2
        yy3d = yy3d.clone() + ny/2
        #zz3d = zz3d.clone() + nz/2

        # Flatten and stack 3D coordinates
        coords3d = torch.stack([xx3d.flatten(), yy3d.flatten(), zz3d.flatten(), torch.ones_like(xx3d).flatten()])
        
        # Get inverse rotation to smear in the opposite direction
        R_inv = compute_affine_matrix(phi, theta, psi, shiftX, shiftY, invertMatrix=True, device=device)
        smeared_coords = R_inv @ coords3d

        # Unflatten smeared coordinates
        smeared_coords = smeared_coords.reshape(4, nx, ny, nz)

        # Check which coordinates are within bounds
        mask = (smeared_coords[0] >= 0) & (smeared_coords[0] < nx) & (smeared_coords[1] >= 0) & (smeared_coords[1] < ny) & (smeared_coords[2] >= 0) & (smeared_coords[2] < nz)

        # Fetch the 2D smeared coordinates
        smeared_coords_2d = smeared_coords[:2].long()

        # Use the smeared 2D coordinates to fetch values from p_torch
        values_to_add = p_torch[smeared_coords_2d[0, mask], smeared_coords_2d[1, mask]]

        # Update the outmap using these values
        outmap[mask] += values_to_add
        
    # Save the final reconstructed volume
    emprove_core.WriteMRC(outmap.cpu().numpy().flatten().tolist(), outputMRC ,nx, ny, nz,1)


def get_rotation_matrix(phi, theta, psi):
    # Convert angles from degrees to radians
    alpha = np.radians(phi)
    beta = np.radians(theta)
    gamma = np.radians(psi)

    ca = np.cos(alpha)
    cb = np.cos(beta)
    cg = np.cos(gamma)
    sa = np.sin(alpha)
    sb = np.sin(beta)
    sg = np.sin(gamma)
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa
    RMatrix = np.array([[cg * cc - sg * sa,            cg * cs + sg * ca,           -cg * sb],
                   [-sg * cc - cg * sa, -sg * cs + cg * ca, sg * sb],
                   [sc, ss,  cb]])
    return RMatrix



from scipy.ndimage import map_coordinates
import math

#from scipy.ndimage import affine_transform
from scipy.fftpack import fft2, ifftshift, ifftn, fftshift

def kaiser_bessel_window(x, alpha):
    """
    x: normalized distance to the center (should be between 0 and 0.5)
    alpha: shape parameter. Higher values yield a wider window.
    """
    # Ensure x is within the valid range
    if x < 0 or x > 0.5:
        return 0
    else:
        # Compute the window value
        value = (1 - (2*x)**2)
        return np.i0(alpha * np.sqrt(value)) / np.i0(alpha)

def create_circle_mask(image_shape, radius):
    h, w = image_shape
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - w/2)**2 + (Y - h/2)**2)
    
    mask = dist_from_center <= radius
    return mask

def insert_image(image2D, fourier2D_mask, outmap, outmapCounter, phi, theta, psi):
    nx, ny, nz = outmap.shape
    
    # Calculate the center of the volume
    center = np.array([nx // 2, ny // 2, nz // 2])
    
    # Calculate the rotation matrix
    rotation_matrix = get_inverse_rotation_matrix(phi, theta, psi)
    # Compute the centered Fourier transform of the 2D image
    fourier_image2D = ifftshift(fft2(image2D))*fourier2D_mask

    # Rotate the 2D image
    #rotated_image = rotate_image(image2D, rotation_matrix, center)
    rotated_image = rotate_image(fourier_image2D, rotation_matrix, center)
    
    # Insert the rotated image into the 3D volume at its center
    for i in range(rotated_image.shape[0]):
        for j in range(rotated_image.shape[1]):
            x, y, z = rotated_image[j, i]
            if 0 <= x < nx and 0 <= y < ny and 0 <= z < nz:
                outmap[int(z), int(y), int(x)] += fourier_image2D[j, i]
                outmapCounter[int(z), int(y), int(x)] += 1 





def get_inverse_rotation_matrix(phi, theta, psi):
    # Convert angles from degrees to radians
    alpha = np.radians(phi)
    beta = np.radians(theta)
    gamma = np.radians(psi)

    ca = np.cos(alpha)
    cb = np.cos(beta)
    cg = np.cos(gamma)
    sa = np.sin(alpha)
    sb = np.sin(beta)
    sg = np.sin(gamma)
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa
    RMatrix = np.array([[cg * cc - sg * sa, cg * cs + sg * ca, -cg * sb],
                        [-sg * cc - cg * sa, -sg * cs + cg * ca, sg * sb],
                        [sc, ss, cb]])
    
    # Return the transpose of the rotation matrix
    return RMatrix.T



def rotate_image(image2D, rotation_matrix, center):
    height, width = image2D.shape
    rotated_image = np.zeros((height, width, 3), dtype=int)
    offset = np.array([height // 2, width // 2, 0])
    for i in range(height):
        for j in range(width):
            point = np.array([i, j, 0]) - offset
            rotated_point = np.dot(rotation_matrix, point) + center
            rotated_image[i, j] = rotated_point
    return rotated_image



#################################
## equalize_images 
emprove_equalize_images = command.add_parser (
    "equalize_images", description="equalize images", help='equalize images'
)
emprove_equalize_images.add_argument("--i", required=True, nargs='+', type=str, help="files with the input MRC maps")
emprove_equalize_images.add_argument("--o_suffix", default='_amplEqualized', type=str, help="suffix to add before the file extension")
emprove_equalize_images.add_argument("--dir", default='./', type=str, help="output directory")


def equalize_images(args):
    fourier_transforms = []
    sizes = []
    sum_amplitudes = None
    count = 0

    # First pass: Compute Fourier transforms and sum their amplitudes
    for input_file in args.i:
        sizeMap = emprove_core.sizeMRC(input_file)
        inputMap = np.array(emprove_core.ReadMRC(input_file))
        inputMap = np.reshape(inputMap, sizeMap)

        spacingMRC = round(emprove_core.spacingMRC(input_file),4)
        #print ("spacing=",spacingMRC)

        ft_inputMap = fftn(inputMap)
        amplitude = np.abs(ft_inputMap)

        if sum_amplitudes is None:
            sum_amplitudes = np.zeros_like(amplitude, dtype=np.complex128)

        sum_amplitudes += amplitude
        fourier_transforms.append(ft_inputMap)
        sizes.append(sizeMap)
        count += 1

    # Calculate average amplitude
    average_amplitude = sum_amplitudes / count
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)


    # Second pass: Replace amplitude and inverse Fourier transform
    for i, ft_map in enumerate(fourier_transforms):
        modified_ft = average_amplitude * (ft_map / np.abs(ft_map))
        modified_map = np.real(ifftn(modified_ft))
        nx, ny, nz = sizes[i]
        input_file = args.i[i]
        base_name, ext = os.path.splitext(os.path.basename(input_file))
        output_filename = f"{base_name}{args.o_suffix}{ext}"
        output_filename = os.path.join(args.dir, output_filename)
        emprove_core.WriteMRC(modified_map.flatten().tolist(), output_filename, nx, ny, nz, spacingMRC)




#################################
## scores_to_csv 
def retrieveScoringTag(columnList):
    def check_sore_tag_correct_format(s):
        # Splitting the string by underscores
        parts = s.split('_')
        # Checking the structure of the string
        if len(parts) == 8 and parts[1] == 'emprove' and parts[2] == 'SCI' and parts[5] == 'scored' and parts[6] == 'selection':
            return True
        return False
    matching_string = None
    for name in columnList:
        if check_sore_tag_correct_format(name):
            matching_string = name
    return matching_string

emprove_scores_to_csv = command.add_parser (
    "scores_to_csv", description="scores_to_csv", help='scores_to_csv'
)
emprove_scores_to_csv.add_argument("--i", required=True, nargs='+', type=str, help="files with the input scores")
emprove_scores_to_csv.add_argument("--csv", default='scoreFile.csv', type=str, help="score_file_to_csv")
emprove_scores_to_csv.add_argument("--o", default='scoreFile.star', type=str, help="scoreFile.star")

def scores_to_csv(args):
    result_df = pd.DataFrame()
    file_class_mapping = {file_name: index + 1 for index, file_name in enumerate(args.i)}
    reference_star=None
    for input_file in args.i:
        if reference_star == None:
            reference_star=input_file
        columns=starHandler.header_columns(input_file)
        tagScore=retrieveScoringTag(columns)
        referenceColumns = starHandler.readColumns(input_file, [tagScore])
        result_df[input_file] = referenceColumns[tagScore]
    result_df['Max_Score_File'] = result_df.idxmax(axis=1)
    result_df['_rlnClassNumber'] = result_df['Max_Score_File'].map(file_class_mapping)
    result_df.to_csv(args.csv)
    #now update the class on the star file
    version=starHandler.infoStarFile(input_file)[2]
    main_section_name="particles"
    if version=="relion_v30":
        main_section_name=""
    starHandler.removeColumnsTagsStartingWith(reference_star, args.o, "_emprove_")
    starHandler.replace_star_columns_from_sections(args.o, args.o, main_section_name,'_rlnClassNumber',result_df)
    




emprove_signal_subtraction_stack = command.add_parser (
    "signal_subtraction_stack", description="signal_subtraction_stack", help='signal_subtraction_stack'
)
emprove_signal_subtraction_stack.add_argument("--star", required=True, type=str, help="star file")
emprove_signal_subtraction_stack.add_argument("--i", required=True, type=str, help="reconstructed map (or fist half map if --i2 is provided)")
emprove_signal_subtraction_stack.add_argument("--i2", required=False, default="None", type=str, help="reconstructed second half maps")
emprove_signal_subtraction_stack.add_argument("--mask", required=False, type=str, help="mrc star file")
emprove_signal_subtraction_stack.add_argument("--o", required=True, type=str, help="output basename")
#emprove_signal_subtraction_stack.add_argument("--angpix", required=True, type=float, help="angpix")

def signal_subtraction_stack(args):
    #multiply mask by maps here
    utils.create_signal_subtraction_stack(args.star, args.i, args.i, args.mask, args.o, saveOriginal=False)







emprove_signal_subtraction_stack = command.add_parser (
    "signal_subtraction_stack", description="signal_subtraction_stack", help='signal_subtraction_stack'
)
emprove_signal_subtraction_stack.add_argument("--star", required=True, type=str, help="star file")
emprove_signal_subtraction_stack.add_argument("--i", required=True, type=str, help="reconstructed map (or fist half map if --i2 is provided)")
emprove_signal_subtraction_stack.add_argument("--i2", required=False, default="None", type=str, help="reconstructed second half maps")
emprove_signal_subtraction_stack.add_argument("--mask", required=False, type=str, help="mrc star file")
emprove_signal_subtraction_stack.add_argument("--o", required=True, type=str, help="output basename")
#emprove_signal_subtraction_stack.add_argument("--angpix", required=True, type=float, help="angpix")

def signal_subtraction_stack(args):
    #multiply mask by maps here
    utils.create_signal_subtraction_stack(args.star, args.i, args.i, args.mask, args.o, saveOriginal=False)





emprove_randomize_halves = command.add_parser (
    "randomize_halves", description="randomize_halves", help='randomize_halves'
)
emprove_randomize_halves.add_argument("--i", required=True, type=str, help="input star file")
emprove_randomize_halves.add_argument("--o", required=True, type=str, help="output starfile")
def randomize_halves_utils(args):
    starHandler.randomize_halves(args.i, args.o)
 



#################################
emprove_extract_particles_from_label_value = command.add_parser (
    "extract_particles_from_label_value", description="extract_particles_from_label_value", help='extract_particles_from_label_value'
)
emprove_extract_particles_from_label_value.add_argument("--i", required=True, type=str, help="input star file")
emprove_extract_particles_from_label_value.add_argument("--label", required=True, type=str, help="label to extract values from")
emprove_extract_particles_from_label_value.add_argument("--value", required=True, type=str, help="value to extract")
emprove_extract_particles_from_label_value.add_argument("--o", required=True, type=str, help="output star file")
def extract_particles_from_label_value(args):
    version=starHandler.infoStarFile(args.i)[2]
    main_section_name="particles"
    if version=="relion_v30":
        main_section_name=""
    starHandler.extract_particles_from_label_from_sections(args.i, args.o, main_section_name, args.label, args.value)


def main(command_line=None):
    args = emprove_parser.parse_args(command_line)
    if args.command == "maskedCrop":
        maskedCrop(args)
    elif args.command == "rotation_average_2d":
        rotation_average_2d(args)
    elif args.command == "rotation_average_stack":
        rotation_average_stack(args)
    elif args.command == "equalize_images":
        equalize_images(args)
    elif args.command == "signal_subtraction_stack":
        signal_subtraction_stack(args)
    elif args.command == "randomize_halves":
        randomize_halves_utils(args)
    elif args.command == "ctf_correct_stack":
        ctf_correct_stack(args)
    elif args.command == "scores_to_csv":
        scores_to_csv(args)
    elif args.command == "extract_particles_from_label_value":
        extract_particles_from_label_value(args)
    else:
        emprove_parser.print_help()


if __name__ == "__main__":
    main()

