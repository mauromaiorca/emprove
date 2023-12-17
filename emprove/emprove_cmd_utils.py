#!/usr/bin/python3


import argparse
import os.path
import numpy as np
import emprove_core
from emprove import starHandler
from emprove import projector_torch
from os import PathLike
from emprove import utils

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





def FBP_reconstruction_cpp2(starFile, outputMRC):
    # Presuming some missing parts of the code based on context
    version=starHandler.infoStarFile(starFile)[2]
    print("star file version=",version)
    #print("sin(180)=",np.sin(math.radians(180)))
    
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
    mask2D=create_circle_mask((nx,ny), nx/2)


    outmapComplex=np.zeros((nx, ny, nz), dtype=complex)
    outmapCounter=np.zeros((nx, ny, nz))
    for index, row in particle_parameters.iterrows():
        tmpLine = row['_rlnImageName']
        atPosition = tmpLine.find('@')
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition+1:]
        phi = (float(row['_rlnAngleRot']))
        theta = (float(row['_rlnAngleTilt']))
        psi = (float(row['_rlnAnglePsi']))
        print("phi=",phi,"  theta=",theta,"     psi=",psi)
        if version=="relion_v31":
            shiftX = float(row['_rlnOriginXAngst'])
            shiftY = float(row['_rlnOriginYAngst'])
        else:
            shiftX = float(row['_rlnOriginX'])
            shiftY = float(row['_rlnOriginY'])
        p=emprove_core.ReadMrcSlice(stackName, imageNo-1)
        image2D=np.reshape(p, (nx, ny))
        insert_image(image2D, mask2D, outmapComplex, outmapCounter, phi, theta, psi)
        #emprove_core.backprojectParticles(p, outmap ,nx, ny, nz,1,phi,theta,psi,shiftX,shiftY)


    #normalize
    outmapComplex[outmapCounter == 0] = 0
    non_zero_indices = outmapCounter != 0
    outmapComplex[non_zero_indices] /= outmapCounter[non_zero_indices]
    outmap_inverse_fourier = np.real(ifftn(fftshift(outmapComplex)))

    #outmap_inverse_fourier_temp = np.real(ifftn(fftshift(outmapComplex)))
    #outmap_inverse_fourier = np.roll(outmap_inverse_fourier_temp, shift=(nx//2, ny//2, nz//2), axis=(0, 1, 2))

    # Save the final reconstructed volume
    emprove_core.WriteMRC(outmap_inverse_fourier.flatten().tolist(), outputMRC ,nx, ny, nz,1)



def FBP_reconstruction_cpp(starFile, outputMRC):
    from scipy.fft import fft2, ifftn, fftshift
    from scipy.interpolate import RegularGridInterpolator
    # Presuming some missing parts of the code based on context
    version=starHandler.infoStarFile(starFile)[2]
    print("star file version=",version)
    #print("sin(180)=",np.sin(math.radians(180)))
    
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
    mask2D=create_circle_mask((nx,ny), nx/2)
    fourier_3d = np.zeros((nz, ny, nx), dtype=complex)


    outmapComplex=np.zeros((nx, ny, nz), dtype=complex)
    outmapCounter=np.zeros((nx, ny, nz))
    counter=0
    for index, row in particle_parameters.iterrows():
        counter+=1
        tmpLine = row['_rlnImageName']
        atPosition = tmpLine.find('@')
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition+1:]
        phi0 = (float(row['_rlnAngleRot']))
        theta0 = (float(row['_rlnAngleTilt']))
        psi0 = (float(row['_rlnAnglePsi']))
        print("phi=",phi0,"  theta=",theta0,"     psi=",psi0)
        if version=="relion_v31":
            shiftX = float(row['_rlnOriginXAngst'])
            shiftY = float(row['_rlnOriginYAngst'])
        else:
            shiftX = float(row['_rlnOriginX'])
            shiftY = float(row['_rlnOriginY'])
        p=emprove_core.ReadMrcSlice(stackName, imageNo-1)
        projection=np.reshape(p, (nx, ny))
        #insert_image(image2D, mask2D, outmapComplex, outmapCounter, phi, theta, psi)
        #emprove_core.backprojectParticles(p, outmap ,nx, ny, nz,1,phi,theta,psi,shiftX,shiftY)
        projection_ft = fft2(projection)
        kz, ky, kx = np.mgrid[:nz, :ny, :nx] - nx // 2
        R = get_inverse_rotation_matrix(phi0, theta0, psi0)
        # Apply a rotation around the z-axis using the angles_psi parameter
        rotated_coords = np.dot(R, np.array([kx.ravel(), ky.ravel(), kz.ravel()]))
        kx_rot, ky_rot, kz_rot = rotated_coords.reshape(3, nz, ny, nx)
        # Create an interpolating function for the 2D Fourier Transform of the projection
        interp_func = RegularGridInterpolator(
            (np.arange(ny), np.arange(nx)),
            projection_ft,
            bounds_error=False,
            fill_value=0
        )

        # Evaluate the interpolating function at the coordinates of the 3D Fourier space
        interpolated_values = interp_func(np.column_stack([ky_rot.ravel(), kx_rot.ravel()]))
        
        # Reshape the interpolated values to match the shape of the 3D Fourier space
        fourier_slice = interpolated_values.reshape((nz, ny, nx))
        
        fourier_3d += fourier_slice

    #normalize
    fourier_3d /= counter
    #outmapComplex[outmapCounter == 0] = 0
    #non_zero_indices = outmapCounter != 0
    #outmapComplex[non_zero_indices] /= outmapCounter[non_zero_indices]
    #outmap_inverse_fourier = np.real(ifftn(fftshift(outmapComplex)))
    reconstructed_object = np.real(fftshift(ifftn(ifftshift(fourier_3d))))

    #outmap_inverse_fourier_temp = np.real(ifftn(fftshift(outmapComplex)))
    #outmap_inverse_fourier = np.roll(outmap_inverse_fourier_temp, shift=(nx//2, ny//2, nz//2), axis=(0, 1, 2))

    # Save the final reconstructed volume
    emprove_core.WriteMRC(reconstructed_object.flatten().tolist(), outputMRC ,nx, ny, nz,1)





emprove_reconstruct_volume = command.add_parser (
    "reconstruct_volume", description="reconstruct a volume, based on the star file", help='reconstruct a volume, based on the star file'
)
emprove_reconstruct_volume.add_argument("--i", required=True, type=str, help="input star file")
emprove_reconstruct_volume.add_argument("--subset", required=False, type=str, default="", help="file with the input mask mrc map")
emprove_reconstruct_volume.add_argument("--mask", required=False, type=str, default="", help="file with the input mask mrc map")
emprove_reconstruct_volume.add_argument("--padding", required=False, type=int, default=2, help="padding value for the mask, default=2")
emprove_reconstruct_volume.add_argument("--o", required=True, type=str, help="output mrc map")
def reconstruct_volume(args):
    #print("ART 3D reconstruction")
    #device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    #ART_3D_reconstruction(args.i, args.o,iterations=1, device=device)
    print("FBP reconstruction")
    #evice = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    #FBP_reconstruction(args.i, args.o,device=device)
    FBP_reconstruction_cpp(args.i, args.o)



emprove_test = command.add_parser (
    "test", description="reconstruct a volume, based on the star file", help='reconstruct a volume, based on the star file'
)
emprove_test.add_argument("--i", required=True, type=str, help="input star file")
emprove_test.add_argument("--o", required=False, type=str, default="", help="output star file")
def test(args):
    print("test")
    import pandas as pd
    df = pd.DataFrame({"_rlnClassNumber": ["1","4","1"]})
    starHandler.update_star_columns_from_sections(args.i, args.o,'particles',df)

#    df = pd.DataFrame({"_rlnVoltage6": ["323"]})
#    starHandler.update_star_columns_from_sections(args.i, args.o,'optics',df)

#    starHandler.delete_star_columns_from_sections(args.i, args.o, 'optics', '_rlnImage')

#    starHandler.delete_star_columns_from_sections(args.i, args.o, 'particles', '_rlnAngle')
#    df = starHandler.merge_star_section(args.i)
#    print (df)






def ParticleVsReprojectionScores_GPU(particlesStarFile: PathLike, scoredParticlesStarFile: PathLike, referenceMap: PathLike, referenceMask: PathLike, angpix, numProcesses = 1, numViews=[50,200,400], doCTF=False):
    mapI=emprove_core.ReadMRC(referenceMap)
    sizeMap=emprove_core.sizeMRC(referenceMap)
    maskI=emprove_core.ReadMRC(referenceMask)



import numpy as np
import scipy.ndimage as ndi



def gaussian_derivative(input_img, sigma, order, axis):
    """
    This function takes an image, performs Gaussian smoothing, and calculates the derivative of the 
    smoothed image. The order of the derivative (1st or 2nd) and the axis along which the derivative is calculated 
    can be specified. The smoothing is done using a Gaussian kernel with a specified standard deviation ('sigma').
    """
    if order not in [0, 1, 2]:
        raise ValueError("Order must be 0, 1 or 2: 0 is blurring, 1 is first derivative, 2 is second derivative")

    # Validate axis
    if axis < 0 or axis >= len(input_img.shape):
        raise ValueError("Invalid axis, 0 is x axis, 1 is y axis, and so on")
    axis=input_img.ndim-axis-1

    # Apply Gaussian smoothing
    smoothed_img = ndi.gaussian_filter(input_img, sigma=sigma, mode='mirror')

    # For 0th order derivative (i.e., just the blur)
    if order == 0:
        return smoothed_img

    # Calculate derivatives
    if order == 1:
        # Calculate first derivative
        # Use np.gradient to calculate the derivative
        # This will automatically handle the differences and keep the shape consistent
        derivatives = np.gradient(smoothed_img, axis=axis)
    elif order == 2:
        # Calculate second derivative
        # Apply the gradient function twice
        first_derivatives = np.gradient(smoothed_img, axis=axis)
        derivatives = np.gradient(first_derivatives, axis=axis)
    return derivatives


import torch
import torch.nn.functional as F

def gaussian_derivative_torch(input_img, sigma, order, axis, device="cpu"):
    """
    This function takes an image, performs Gaussian smoothing, and calculates the derivative of the 
    smoothed image. The order of the derivative (1st or 2nd) and the axis along which the derivative is calculated 
    can be specified. The smoothing is done using a Gaussian kernel with a specified standard deviation ('sigma').
    """
    if order not in [0, 1, 2]:
        raise ValueError("Order must be 0, 1 or 2: 0 is blurring, 1 is first derivative, 2 is second derivative")

    # Validate axis
    if axis < 0 or axis >= len(input_img.shape):
        raise ValueError("Invalid axis, 0 is x axis, 1 is y axis, and so on")

    # PyTorch's gaussian_filter equivalent
    smoothed_img = F.gaussian_blur(input_img, kernel_size=int(3*sigma), sigma=sigma)

    # For 0th order derivative (i.e., just the blur)
    if order == 0:
        return smoothed_img

    # Calculate derivatives
    if order == 1:
        # Calculate first derivative
        derivatives = torch.gradient(smoothed_img, dim=axis)
    elif order == 2:
        # Calculate second derivative
        first_derivatives = torch.gradient(smoothed_img, dim=axis)
        derivatives = torch.gradient(first_derivatives, dim=axis)

    return derivatives


def cross_correlation(image1, image2, mask=None):
    """
    This function calculates the normalized cross-correlation between two n-dimensional images, 
    given a mask. The mask is used to ignore the corresponding points in the images from the calculation of 
    cross-correlation. The function returns a single correlation value.
    """

    assert image1.shape == image2.shape, "Images must be the same size"
    result = 0
    norm1 = 0
    norm2 = 0

    if mask is not None:
        mask = np.array(mask, dtype=bool)

    sum1 = 0
    sum2 = 0
    count = 0

    for index, value in np.ndenumerate(image1):
        if mask is None or mask[index]:
            sum1 += image1[index]
            sum2 += image2[index]
            count += 1

    mean1 = sum1 / count
    mean2 = sum2 / count

    # Subtract means from images and compute the cross-correlation, 
    # considering the mask if provided
    for index, value in np.ndenumerate(image1):
        if mask is None or mask[index]:
            val1 = image1[index] - mean1
            val2 = image2[index] - mean2
            result += val1 * val2
            norm1 += val1**2
            norm2 += val2**2

    # Normalize result
    result /= np.sqrt(norm1 * norm2)

    return max(0,result)


def cross_correlation_torch(image1, image2, mask, device="cpu"):
    """
    This function calculates the normalized cross-correlation between two n-dimensional images, 
    given a mask. The mask is used to ignore the corresponding points in the images from the calculation of 
    cross-correlation. The function returns a single correlation value.
    """
    assert image1.shape == image2.shape, "Images must be the same size"
    sum1 = torch.sum(image1[mask]) if mask is not None else torch.sum(image1)
    sum2 = torch.sum(image2[mask]) if mask is not None else torch.sum(image2)
    count = torch.sum(mask) if mask is not None else torch.numel(image1)
    mean1 = sum1 / count
    mean2 = sum2 / count
    val1 = image1 - mean1
    val2 = image2 - mean2
    result = torch.sum(val1 * val2)
    norm1 = torch.sum(val1**2)
    norm2 = torch.sum(val2**2)
    result /= torch.sqrt(norm1 * norm2)
    return torch.clamp(result, min=0)



def SCI(I, RI, sigma=1.0, MaskImage=None):
    """
    This function calculates the Structural Cross-Correlation Index (SCI) between two n-dimensional images. 
    The SCI is a measure of the similarity between two images. The process involves performing a Fourier transform 
    on the images, calculating an average amplitude spectrum, applying that to the original images, and then 
    calculating cross-correlation of the original and derivative images.
    """
    assert I.shape == RI.shape, "Images must be the same size"
    I_fft = np.fft.fftn(I)
    RI_fft = np.fft.fftn(RI)
    I_abs_fft = np.abs(I_fft) + 1e-7
    RI_abs_fft = np.abs(RI_fft) + 1e-7
    ampAvg = 0.5 * (I_abs_fft + RI_abs_fft)
    I_out = np.fft.ifftn(I_fft * (ampAvg / I_abs_fft)).real
    RI_out = np.fft.ifftn(RI_fft * (ampAvg / RI_abs_fft)).real

    ndim = I_out.ndim
    scoreSCI_final = cross_correlation(RI_out, I_out, MaskImage)
    for axis in range(ndim):
        I_out_d1 = gaussian_derivative(I_out, sigma, 1, axis)
        RI_out_d1 = gaussian_derivative(RI_out, sigma, 1, axis)
        I_out_d2 = gaussian_derivative(I_out, sigma, 2, axis)
        RI_out_d2 = gaussian_derivative(RI_out, sigma, 2, axis)
        scoreSCI_final *= cross_correlation(RI_out_d1, I_out_d1, MaskImage) 
        scoreSCI_final *= cross_correlation(RI_out_d2, I_out_d2, MaskImage)
    return scoreSCI_final

def SCI_torch(I, RI, MaskImage, sigma=1.0, device="cpu"):
    """
    This function calculates the Structural Cross-Correlation Index (SCI) between two n-dimensional images. 
    The SCI is a measure of the similarity between two images. The process involves performing a Fourier transform 
    on the images, calculating an average amplitude spectrum, applying that to the original images, and then 
    calculating cross-correlation of the original and derivative images.
    """
    assert I.shape == RI.shape, "Images must be the same size"
    I = I.to(device)
    RI = RI.to(device)
    #if MaskImage is not None:
    #    MaskImage = MaskImage.to(device)

    I_fft = torch.fft.fftn(I)
    RI_fft = torch.fft.fftn(RI)
    I_abs_fft = torch.abs(I_fft) + 1e-7
    RI_abs_fft = torch.abs(RI_fft) + 1e-7
    ampAvg = 0.5 * (I_abs_fft + RI_abs_fft)
    I_out = torch.fft.ifftn(I_fft * (ampAvg / I_abs_fft)).real
    RI_out = torch.fft.ifftn(RI_fft * (ampAvg / RI_abs_fft)).real

    ndim = I_out.ndim
    scoreSCI_final = cross_correlation_torch(RI_out, I_out, MaskImage, device)
    for axis in range(ndim):
        I_out_d1 = gaussian_derivative_torch(I_out, sigma, 1, axis, device)
        RI_out_d1 = gaussian_derivative_torch(RI_out, sigma, 1, axis, device)
        I_out_d2 = gaussian_derivative_torch(I_out, sigma, 2, axis, device)
        RI_out_d2 = gaussian_derivative_torch(RI_out, sigma, 2, axis, device)
        scoreSCI_final *= cross_correlation_torch(RI_out_d1, I_out_d1, MaskImage, device)
        scoreSCI_final *= cross_correlation_torch(RI_out_d2, I_out_d2, MaskImage, device)
    return scoreSCI_final



emprove_assess_particles_gpu = command.add_parser (
    "assess_particles_gpu", description="assess_particles_gpu", help='assess_particles_gpu'
)
emprove_assess_particles_gpu.add_argument("--star", required=False, type=str, help="input star file")
emprove_assess_particles_gpu.add_argument("--o", required=False, type=str, default="", help="output star file")
emprove_assess_particles_gpu.add_argument("--map", required=True, type=str, default="", help="referenceMap map for comparison")
emprove_assess_particles_gpu.add_argument("--mask", required=True, type=str, default="", help="mask")
emprove_assess_particles_gpu.add_argument("--angpix", required=False, type=str, default="", help="output star file")

def assess_particles_gpu(args):
    print("emprove_assess_particles_gpu")
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")  # use GPU 0 if available, else CPU

    angpix=utils.get_MRC_map_pixel_spacing(args.map)
    print ("angpix=",angpix[0])

    #read the starFile
    particlesStarFile=args.star
    version=starHandler.infoStarFile(particlesStarFile)[2]
    print("star file version:",version)

    if version=="relion_v31":
        parameters=starHandler.merge_star_section(particlesStarFile, nameMainSection="particles", nameIndexedSection="optics", indexingTag="_rlnOpticsGroup")
        #print ("parameters before")
        #print (parameters)
        tmp_angpix=float(angpix[0])
        if '_rlnImagePixelSize' in parameters.columns:
            parameters['_rlnImagePixelSize'] = parameters['_rlnImagePixelSize'].astype(float)
            tmp_angpix=parameters['_rlnImagePixelSize']
        else:
            parameters['_rlnImagePixelSize'] = tmp_angpix
        if '_rlnOriginX' not in parameters.columns and '_rlnOriginXAngst' in parameters.columns:
            parameters['_rlnOriginX'] = parameters['_rlnOriginXAngst'].astype(float)/tmp_angpix
            parameters=parameters.drop(['_rlnOriginXAngst'],axis=1)
        if '_rlnOriginY' not in parameters.columns and '_rlnOriginYAngst' in parameters.columns:
            parameters['_rlnOriginY']=parameters['_rlnOriginYAngst'].astype(float)/tmp_angpix
            parameters=parameters.drop(['_rlnOriginYAngst'],axis=1)
        #print ("parameters after")
        #print (parameters)
    else:
        parameters=starHandler.readStar(particlesStarFile)
        parameters['_rlnOriginX']=parameters['_rlnOriginX'].astype(float)
        parameters['_rlnOriginY']=parameters['_rlnOriginY'].astype(float)
        if '_rlnImagePixelSize' not in parameters.columns:
            parameters['_rlnImagePixelSize'] = float(angpix[0])
    if '_rlnPhaseShift' not in parameters.columns:
        parameters['_rlnPhaseShift'] = 0
    numParticles=len(parameters['_rlnImageName'])

    doCtf=False
    #############
    #do the score_block_particles

    #load map to memory
    sizeMap=emprove_core.sizeMRC(args.map)
    inputMap=np.array(emprove_core.ReadMRC(args.map))
    inputMap=np.reshape(inputMap,(sizeMap[2],sizeMap[1],sizeMap[0]))
    inputMapTorch = torch.tensor(inputMap, dtype=torch.float32)

    #load mask to memory
    sizeMask=emprove_core.sizeMRC(args.mask)
    inputMask=(np.array(emprove_core.ReadMRC(args.mask))> 0.2)
    inputMask=np.reshape(inputMask,(sizeMask[2],sizeMask[1],sizeMask[0]))
    inputMaskTorch = torch.tensor(inputMask, dtype=torch.bool)

    emprove_core.WriteEmptyMRC("reprojection_stack.mrc",sizeMap[0],sizeMap[1],numParticles)
    emprove_core.WriteEmptyMRC("reprojection_stack_torch.mrc",sizeMap[0],sizeMap[1],numParticles)
    import time
    import assessParticles_GPU

    #
    for index, row in parameters.iterrows():
        print("particle ", index)

        #parameters Position
        phi = torch.tensor(np.radians(float(row['_rlnAngleRot'])), dtype=torch.float32).to(device)
        theta = torch.tensor(np.radians(float(row['_rlnAngleTilt'])), dtype=torch.float32).to(device)
        psi =  torch.tensor(np.radians(float(row['_rlnAnglePsi'])), dtype=torch.float32).to(device)



        #extract the slice
        imageNo=int(row['_rlnImageName'].split('@')[0]) - 1
        stackName=row['_rlnImageName'].split('@')[1]
        I=emprove_core.ReadMrcSlice(stackName,imageNo)
        
        #parameters Ctf
        if doCtf:
            Voltage=float(row['_rlnVoltage'])
            DefocusU=float(row['_rlnDefocusU']) 
            DefocusV=float(row['_rlnDefocusV'])
            DefocusAngle=float(row['_rlnDefocusAngle']) 
            SphericalAberration=float(row['_rlnSphericalAberration'])
            CtfBfactor=float(row['_rlnCtfBfactor'])
            PhaseShift=float(row['_rlnPhaseShift'])
            AmplitudeContrast=float(row['_rlnAmplitudeContrast'])
            ImagePixelSize=float(row['_rlnImagePixelSize'])
            start_time = time.time()
            I_torch=assessParticles_GPU.transformCtfImage_torch(I,sizeMap[0],sizeMap[1], ImagePixelSize, Voltage, DefocusU, DefocusV, DefocusAngle, SphericalAberration, CtfBfactor, PhaseShift, AmplitudeContrast, device=device)
            end_time = time.time()
            print(f"CUDA executed in : {end_time - start_time} seconds")
            start_time = time.time()
            I=assessParticles_GPU.transformCtfImage(I,sizeMap[0],sizeMap[1], ImagePixelSize, Voltage, DefocusU, DefocusV, DefocusAngle, SphericalAberration, CtfBfactor, PhaseShift, AmplitudeContrast)
            end_time = time.time()
            print(f"c++ executed in : {end_time - start_time} seconds")
        else:
            I_torch=torch.tensor(I, device=device).reshape(sizeMap[0],sizeMap[1])


        #project the map
        start_time = time.time()
        projection = projector_torch.forward_projection(inputMapTorch, phi, theta, psi, float(row['_rlnOriginX']), float(row['_rlnOriginY']), device="cpu")
        projectionOut = projection.cpu().numpy()
        end_time = time.time()
        print(f"torch CPU projecting executed in : {end_time - start_time} seconds")

        start_time = time.time()
        projection = projector_torch.forward_projection(inputMapTorch, phi, theta, psi, float(row['_rlnOriginX']), float(row['_rlnOriginY']), device="cuda:0")
        projectionOut = projection.cpu().numpy()
        end_time = time.time()
        print(f"torch CUDA GPU:0 projecting executed in : {end_time - start_time} seconds")

        emprove_core.ReplaceMrcSlice(projectionOut.flatten().tolist(),"reprojection_stack_torch.mrc",sizeMap[0],sizeMap[1],index)

        start_time = time.time()
        RI=emprove_core.projectMap(inputMap.flatten().tolist(),sizeMap[0],sizeMap[1],sizeMap[2],float(row['_rlnAngleRot']),float(row['_rlnAngleTilt']),float(row['_rlnAnglePsi']),float(row['_rlnOriginX']), float(row['_rlnOriginY']),0)
        end_time = time.time()
        print(f"c++ projecting executed in : {end_time - start_time} seconds")
        

        MI=emprove_core.projectMask(inputMask.flatten().tolist(),sizeMap[0],sizeMap[1],sizeMap[2],float(row['_rlnAngleRot']),float(row['_rlnAngleTilt']),float(row['_rlnAnglePsi']),float(row['_rlnOriginX']), float(row['_rlnOriginY']),0,0.5)
        emprove_core.ReplaceMrcSlice(RI,"reprojection_stack.mrc",sizeMap[0],sizeMap[1],index)

        #comparison:

        #c++
        start_time = time.time()
        I_fft=np.fft.fftn(np.reshape(I,[sizeMap[1],sizeMap[0]]))
        RI_fft=np.fft.fftn(np.reshape(RI,[sizeMap[1],sizeMap[0]]))
        I_abs_fft=np.abs(I_fft)+0.0000001
        RI_abs_fft=np.abs(RI_fft)+0.0000001
        ampAvg=0.5*(I_abs_fft+RI_abs_fft)
        I_out=np.fft.ifftn(I_fft*(ampAvg/I_abs_fft)).real.flatten().tolist()
        RI_out=np.fft.ifftn(RI_fft*(ampAvg/RI_abs_fft)).real.flatten().tolist()
        scoreSCI=emprove_core.MaskedImageComparison(RI_out, I_out, MI, sizeMap[0], sizeMap[1], 1,"SCI","unprocessed","1")
        print ("c++ SCI score=",scoreSCI)
        end_time = time.time()
        print(f"c++ SCI executed in : {end_time - start_time} seconds")

        #numpy
        start_time = time.time()
        scoreSCI_NP=SCI(np.reshape(I,[sizeMap[1],sizeMap[0]]), np.reshape(RI,[sizeMap[1],sizeMap[0]]), sigma=1.0, MaskImage=np.reshape(MI,[sizeMap[1],sizeMap[0]]))
        end_time = time.time()
        print ("scoreSCI_NP score=",scoreSCI_NP)
        print(f"scoreSCI_NP executed in : {end_time - start_time} seconds")

        #torch
        #start_time = time.time()
        #M_torch=torch.tensor(np.reshape(MI,[sizeMap[1],sizeMap[0]]), dtype=torch.bool).to(device)
        #scoreSCI_torch=SCI_torch(I_torch, projection, M_torch, sigma=1.0, device=device)
        #end_time = time.time()
        #print ("scoreSCI_torch score=",scoreSCI_torch)
        #print(f"scoreSCI_torch executed in : {end_time - start_time} seconds")



        #compare it


    return
    doCTF=True
    if doCTF:
        columns=starHandler.header_columns(particlesStarFile)
        print (len(coordinates))
        if not '_rlnPhaseShift' in columns:
            PhaseShift=pd.DataFrame(np.zeros(len(coordinates)))
        else:
            PhaseShift = starHandler.readColumns(particlesStarFile, ['_rlnPhaseShift'])

        #print (columns)

        print ("doing CTF...")
        if version=="relion_v31":
            print('READ parameters')
            parametersFULL = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnOpticsGroup','_rlnCtfBfactor'])
            print('GOT the parameters')
            idx=[x for x in range(0, len(parametersFULL))]
            parametersFULL['idx']=idx
            parametersDataOptics = starHandler.dataOptics(particlesStarFile)[['_rlnImagePixelSize','_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration','_rlnOpticsGroup']]
            ctfParameters =  pd.merge(parametersFULL, parametersDataOptics,  on=['_rlnOpticsGroup']).sort_values(['idx'])
            ctfParameters=ctfParameters.drop(['_rlnOpticsGroup'],axis=1).reindex()
            ctfParameters=ctfParameters.set_index('idx')
            ctfParameters.rename(columns={'_rlnImagePixelSize':'_rlnDetectorPixelSize'},inplace=True)
        else:
            (ctfParameters) = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnDetectorPixelSize','_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration','_rlnCtfBfactor'])
        ctfParameters['_rlnPhaseShift']=PhaseShift
    numParticles=len(coordinates['_rlnImageName'])
    print('num particles=',numParticles,'  num Processes=',numProcesses)

    if numProcesses > multiprocessing.cpu_count():
        numProcesses = multiprocessing.cpu_count()
    elif numProcesses < 1:
        numProcesses = 1
    if numProcesses > numParticles:
        numProcesses = numParticles

    blockSize=(np.floor(numParticles/numProcesses))
    idxMatrix=np.zeros([numProcesses,2])
    for mm in range (0, numProcesses):
        idxMatrix[mm][0]=(mm*blockSize)
        idxMatrix[mm][1]=((mm+1)*(blockSize))
        if mm == numProcesses-1:
            idxMatrix[mm][1]=int(numParticles)


    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs=[]
    for ii in range(0, len(idxMatrix) ):
        #print ('**************\n***********\ndebug=',idxMatrix[ii][0]),'   ',int(idxMatrix[ii][1])
        subsetCtfParameters=ctfParameters[int(idxMatrix[ii][0]):int(idxMatrix[ii][1])]
        p=multiprocessing.Process(target=scoreBlockParticles, args=(idxMatrix[ii][0],idxMatrix[ii][1],listScoresTags, mapI, maskI, sizeMap, angpix, coordinates, subsetCtfParameters, ii,return_dict,))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()

    scoresOut=[]
    for ii in range(0, len(idxMatrix) ):
        scoresOut=scoresOut+return_dict.values()[ii].tolist()
    scoresOut = np.array(scoresOut).reshape( [numParticles,len(listScoresTags)+1] )
    listScoreOut = np.array(sorted(scoresOut, key=operator.itemgetter(0), reverse=False)).reshape( [numParticles,len(listScoresTags)+1] )[:,1:].tolist()
    df_full_scores = pd.DataFrame(data=listScoreOut,columns=listScoresTags)
    starHandler.removeColumnsTagsStartingWith(particlesStarFile, scoredParticlesStarFile, "_emprove_")
    starHandler.addDataframeColumns(scoredParticlesStarFile, scoredParticlesStarFile, listScoresTags, df_full_scores)



def main(command_line=None):
    args = emprove_parser.parse_args(command_line)
    if args.command == "maskedCrop":
        maskedCrop(args)
    elif args.command == "rotation_average_2d":
        rotation_average_2d(args)
    elif args.command == "rotation_average_stack":
        rotation_average_stack(args)
    elif args.command == "project_volume":
        project_volume(args)
    elif args.command == "reconstruct_volume":
        reconstruct_volume(args)
    elif args.command == "assess_particles_gpu":
        assess_particles_gpu(args)
    elif args.command == "test":
        test(args)
    else:
        emprove_parser.print_help()


if __name__ == "__main__":
    main()

