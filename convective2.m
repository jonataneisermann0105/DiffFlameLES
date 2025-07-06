function [convu, convv, convw, convZ] = convective2(ni, nj, nk, x, y, z, ra, rua, rva, rwa, rZa)
    
    % Inputs:
    %   ni   - Number of grid points in the x-direction
    %   nj   - Number of grid points in the y-direction
    %   nk   - Number of grid points in the z-direction
    %   x    - 1D array of spatial coordinates in the x-direction
    %   y    - 1D array of spatial coordinates in the y-direction
    %   z    - 1D array of spatial coordinates in the z-direction
    %   ra   - 3D array of fluid density:         ra  = ρ
    %   rua  - 3D array of momentum in x:         rua = ρu
    %   rva  - 3D array of momentum in y:         rva = ρv
    %   rwa  - 3D array of momentum in z:         rwa = ρw
    %   rZa  - 3D array of a conserved scalar:    rZa = ρZ

    % Outputs:
    %   convu - Convective term of the x-momentum equation
    %   convv - Convective term of the y-momentum equation
    %   convw - Convective term of the z-momentum equation
    %   convZ - Convective term of the mixture fraction equation

    % Author: Jonatan Ismael Eisermann
    % Date: July 6, 2025.

    convu = zeros(ni,nj,nk);
    convv = zeros(ni,nj,nk);
    convw = zeros(ni,nj,nk);
    convZ = zeros(ni,nj,nk);

    for k = 2 : nk-1
        for j = 2 : nj-1
            for i = 2 : ni-1
                dx          = 0.5 * (x(i+1) - x(i-1));
                dy          = 0.5 * (y(j+1) - y(j-1));
                dz          = 0.5 * (z(k+1) - z(k-1));
                dx2         = dx * dx;
                dy2         = dy * dy;
                dz2         = dz * dz;

                % Convective terms
                druudx = (rua(i+1,j,k) * rua(i+1,j,k) / ra(i+1,j,k) - rua(i-1,j,k) * rua(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);
                druvdx = (rua(i+1,j,k) * rva(i+1,j,k) / ra(i+1,j,k) - rua(i-1,j,k) * rva(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);
                druwdx = (rua(i+1,j,k) * rwa(i+1,j,k) / ra(i+1,j,k) - rua(i-1,j,k) * rwa(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);
                druZdx = (rua(i+1,j,k) * rZa(i+1,j,k) / ra(i+1,j,k) - rua(i-1,j,k) * rZa(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);

                druvdy = (rua(i,j+1,k) * rva(i,j+1,k) / ra(i,j+1,k) - rua(i,j-1,k) * rva(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);
                drvvdy = (rva(i,j+1,k) * rva(i,j+1,k) / ra(i,j+1,k) - rva(i,j-1,k) * rva(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);
                drvwdy = (rva(i,j+1,k) * rwa(i,j+1,k) / ra(i,j+1,k) - rva(i,j-1,k) * rwa(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);
                drvZdy = (rva(i,j+1,k) * rZa(i,j+1,k) / ra(i,j+1,k) - rva(i,j-1,k) * rZa(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);

                druwdz = (rua(i,j,k+1) * rwa(i,j,k+1) / ra(i,j,k+1) - rua(i,j,k-1) * rwa(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);
                drvwdz = (rva(i,j,k+1) * rwa(i,j,k+1) / ra(i,j,k+1) - rva(i,j,k-1) * rwa(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);
                drwwdz = (rwa(i,j,k+1) * rwa(i,j,k+1) / ra(i,j,k+1) - rwa(i,j,k-1) * rwa(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);
                drwZdz = (rwa(i,j,k+1) * rZa(i,j,k+1) / ra(i,j,k+1) - rwa(i,j,k-1) * rZa(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);

                convu(i,j,k) = druudx + druvdy + druwdz;
                convv(i,j,k) = druvdx + drvvdy + drvwdz;
                convw(i,j,k) = druwdx + drvwdy + drwwdz;
                convZ(i,j,k) = druZdx + drvZdy + drwZdz;
               
            end
        end
    end

    for k = 1 : nk
        for j = 1 : nj
            convu(ni,j,k) = convu(ni-1,j,k);
            convv(ni,j,k) = convv(ni-1,j,k);
            convw(ni,j,k) = convw(ni-1,j,k);
            convZ(ni,j,k) = convZ(ni-1,j,k);
        end
    end

end