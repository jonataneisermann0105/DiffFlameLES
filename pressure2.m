function [dpdx, dpdy, dpdz] = pressure2(ni, nj, nk, x, y, z, pa)

    % Inputs:
    %   ni   - Number of grid points in the x-direction
    %   nj   - Number of grid points in the y-direction
    %   nk   - Number of grid points in the z-direction
    %   x    - 1D array of spatial coordinates in the x-direction
    %   y    - 1D array of spatial coordinates in the y-direction
    %   z    - 1D array of spatial coordinates in the z-direction
    %   pa   - 3D array of pressure field

    % Outputs:
    %   dpdx - Pressure gradient in the x-direction: ∂p/∂x
    %   dpdy - Pressure gradient in the y-direction: ∂p/∂y
    %   dpdz - Pressure gradient in the z-direction: ∂p/∂z

    % Author: Jonatan Ismael Eisermann
    % Date: July 6, 2025.
    
    dpdx = zeros(ni, nj, nk);
    dpdy = zeros(ni, nj, nk);
    dpdz = zeros(ni, nj, nk);

    for k = 2 : nk-1
        for j = 2 : nj-1
            for i = 2 : ni-1
                dx  = (x(i+1) - x(i-1)) / 2;
                dy  = (y(j+1) - y(j-1)) / 2;
                dz  = (z(k+1) - z(k-1)) / 2;

                dpdx(i,j,k) = (pa(i+1,j,k) - pa(i-1,j,k)) / (2 * dx);
                dpdy(i,j,k) = (pa(i,j+1,k) - pa(i,j-1,k)) / (2 * dy);
                dpdz(i,j,k) = (pa(i,j,k+1) - pa(i,j,k-1)) / (2 * dz);
            end
        end
    end

end