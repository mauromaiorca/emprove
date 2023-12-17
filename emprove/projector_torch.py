import torch

def compute_affine_matrix(phi, theta, psi, x_shift=0, y_shift=0, invertMatrix=True, device='cpu'):
    # Angles already converted from degrees to radians
    #phi = torch.tensor(np.radians(phi), dtype=torch.float32).to(device)
    #theta = torch.tensor(np.radians(theta), dtype=torch.float32).to(device)
    #psi = torch.tensor(np.radians(psi), dtype=torch.float32).to(device)

    # Compute trig values
    ca = torch.cos(phi)
    cb = torch.cos(theta)
    cg = torch.cos(psi)
    sa = torch.sin(phi)
    sb = torch.sin(theta)
    sg = torch.sin(psi)
    
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa

    # Construct the matrix using rotation and translation components
    matrix = torch.tensor([
        [cg * cc - sg * sa,    cg * cs + sg * ca,   -cg * sb,  x_shift],
        [-sg * cc - cg * sa,  -sg * cs + cg * ca,   sg * sb,  y_shift],
        [sc,                  ss,                   cb,       0],
        [0,                   0,                    0,        1]
    ], device=device)

    if invertMatrix:
        matrix = torch.inverse(matrix)
    print ("matrix=",matrix)
    return matrix


def forward_projection(volume, phi, theta, psi, x_shift=0, y_shift=0, device='cpu'):
    normalized_x_shift = -2 * x_shift / volume.shape[-1]  # assuming width is the last dimension
    normalized_y_shift = -2 * y_shift / volume.shape[-2]  # assuming height is the second to last dimension
    affine_matrix = compute_affine_matrix(phi, theta, psi,normalized_x_shift, normalized_y_shift, device=device)
    
    # Remove the last row from the 4x4 matrix to have a 3x4 matrix
    affine_matrix = affine_matrix[:-1].clone()  # use clone() to avoid in-place modifications
    if torch.isnan(affine_matrix).any():
        print(f"NaN detected in 'affine_matrix' at iteration {iteration}, index {index}")

    
    # The shape of the volume should be [batch_size, channels, depth, height, width]
    # For our case, we can assume batch_size=1 and channels=1
    volume = volume.unsqueeze(0).unsqueeze(0).to(device).float()

    # Generate the sampling grid
    grid = torch.nn.functional.affine_grid(affine_matrix.unsqueeze(0), volume.size(), align_corners=True).to(device).float()
    
    # Warp the volume based on the grid
    transformed_volume = torch.nn.functional.grid_sample(volume, grid, mode='bilinear', padding_mode='zeros', align_corners=True)
    
    # Remove added batch and channel dimensions
    transformed_volume = transformed_volume.squeeze(0).squeeze(0)

    # Given the parallel beam geometry, sum along the y-axis
    projection = torch.sum(transformed_volume, axis=0)
    return projection

