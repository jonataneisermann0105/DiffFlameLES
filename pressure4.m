function [dpdx, dpdy, dpdz] = pressure4(ni, nj, nk, x, y, z, pa)

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

    % Pressure gradient
    for k = 2 : nk-1
        for j = 2 : nj-1
            for i = 2 : ni-1
                dx  = (x(i+1) - x(i-1)) / 2;
                dy  = (y(j+1) - y(j-1)) / 2;
                dz  = (z(k+1) - z(k-1)) / 2;

                iw = i-2;
                ie = i+2;
                js = j-2;
                jn = j+2;
                kb = k-2;
                kf = k+2;

                if i == 2
                    iw = i-1;
                end
                
                if i == (ni-1)
                    ie = ni;
                end

                if j == 2
                    js = j-1;
                end

                if j == (nj-1)
                    jn = nj;
                end

                if k == 2
                    kb = k-1;
                end

                if k == (nk-1)
                    kf = nk;
                end

                dpdx(i,j,k) = (- pa(ie,j,k) + 8 * pa(i+1,j,k) - 8 * pa(i-1,j,k) + pa(iw,j,k)) / (12 * dx);
                dpdy(i,j,k) = (- pa(i,jn,k) + 8 * pa(i,j+1,k) - 8 * pa(i,j-1,k) + pa(i,js,k)) / (12 * dy);
                dpdz(i,j,k) = (- pa(i,j,kf) + 8 * pa(i,j,k+1) - 8 * pa(i,j,k-1) + pa(i,j,kb)) / (12 * dz);
            end
        end
    end

end