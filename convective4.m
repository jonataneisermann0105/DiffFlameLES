function [convu, convv, convw, convZ] = convective4(ni, nj, nk, x, y, z, ra, rua, rva, rwa, rZa)

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
                dx          = (x(i+1) - x(i-1)) / 2;
                dy          = (y(j+1) - y(j-1)) / 2;
                dz          = (z(k+1) - z(k-1)) / 2;

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

                % Convective terms
                druudx = (- rua(ie,j,k) * rua(ie,j,k) / ra(ie,j,k) + 8 * rua(i+1,j,k) * rua(i+1,j,k) / ra(i+1,j,k) - 8 * rua(i-1,j,k) * rua(i-1,j,k) / ra(i-1,j,k) + rua(iw,j,k) * rua(iw,j,k) / ra(iw,j,k)) / (12 * dx);
                druvdx = (- rua(ie,j,k) * rva(ie,j,k) / ra(ie,j,k) + 8 * rua(i+1,j,k) * rva(i+1,j,k) / ra(i+1,j,k) - 8 * rua(i-1,j,k) * rva(i-1,j,k) / ra(i-1,j,k) + rua(iw,j,k) * rva(iw,j,k) / ra(iw,j,k)) / (12 * dx);
                druwdx = (- rua(ie,j,k) * rwa(ie,j,k) / ra(ie,j,k) + 8 * rua(i+1,j,k) * rwa(i+1,j,k) / ra(i+1,j,k) - 8 * rua(i-1,j,k) * rwa(i-1,j,k) / ra(i-1,j,k) + rua(iw,j,k) * rwa(iw,j,k) / ra(iw,j,k)) / (12 * dx);
                druZdx = (- rua(ie,j,k) * rZa(ie,j,k) / ra(ie,j,k) + 8 * rua(i+1,j,k) * rZa(i+1,j,k) / ra(i+1,j,k) - 8 * rua(i-1,j,k) * rZa(i-1,j,k) / ra(i-1,j,k) + rua(iw,j,k) * rZa(iw,j,k) / ra(iw,j,k)) / (12 * dx);

                druvdy = (- rua(i,jn,k) * rva(i,jn,k) / ra(i,jn,k) + 8 * rua(i,j+1,k) * rva(i,j+1,k) / ra(i,j+1,k) - 8 * rua(i,j-1,k) * rva(i,j-1,k) / ra(i,j-1,k) + rua(i,js,k) * rva(i,js,k) / ra(i,js,k)) / (12 * dy);
                drvvdy = (- rva(i,jn,k) * rva(i,jn,k) / ra(i,jn,k) + 8 * rva(i,j+1,k) * rva(i,j+1,k) / ra(i,j+1,k) - 8 * rva(i,j-1,k) * rva(i,j-1,k) / ra(i,j-1,k) + rva(i,js,k) * rva(i,js,k) / ra(i,js,k)) / (12 * dy);
                drvwdy = (- rva(i,jn,k) * rwa(i,jn,k) / ra(i,jn,k) + 8 * rva(i,j+1,k) * rwa(i,j+1,k) / ra(i,j+1,k) - 8 * rva(i,j-1,k) * rwa(i,j-1,k) / ra(i,j-1,k) + rva(i,js,k) * rwa(i,js,k) / ra(i,js,k)) / (12 * dy);
                drvZdy = (- rva(i,jn,k) * rZa(i,jn,k) / ra(i,jn,k) + 8 * rva(i,j+1,k) * rZa(i,j+1,k) / ra(i,j+1,k) - 8 * rva(i,j-1,k) * rZa(i,j-1,k) / ra(i,j-1,k) + rva(i,js,k) * rZa(i,js,k) / ra(i,js,k)) / (12 * dy);

                druwdz = (- rua(i,j,kf) * rwa(i,j,kf) / ra(i,j,kf) + 8 * rua(i,j,k+1) * rwa(i,j,k+1) / ra(i,j,k+1) - 8 * rua(i,j,k-1) * rwa(i,j,k-1) / ra(i,j,k-1) + rua(i,j,kb) * rwa(i,j,kb) / ra(i,j,kb)) / (12 * dz);
                drvwdz = (- rva(i,j,kf) * rwa(i,j,kf) / ra(i,j,kf) + 8 * rva(i,j,k+1) * rwa(i,j,k+1) / ra(i,j,k+1) - 8 * rva(i,j,k-1) * rwa(i,j,k-1) / ra(i,j,k-1) + rva(i,j,kb) * rwa(i,j,kb) / ra(i,j,kb)) / (12 * dz);
                drwwdz = (- rwa(i,j,kf) * rwa(i,j,kf) / ra(i,j,kf) + 8 * rwa(i,j,k+1) * rwa(i,j,k+1) / ra(i,j,k+1) - 8 * rwa(i,j,k-1) * rwa(i,j,k-1) / ra(i,j,k-1) + rwa(i,j,kb) * rwa(i,j,kb) / ra(i,j,kb)) / (12 * dz);
                drwZdz = (- rwa(i,j,kf) * rZa(i,j,kf) / ra(i,j,kf) + 8 * rwa(i,j,k+1) * rZa(i,j,k+1) / ra(i,j,k+1) - 8 * rwa(i,j,k-1) * rZa(i,j,k-1) / ra(i,j,k-1) + rwa(i,j,kb) * rZa(i,j,kb) / ra(i,j,kb)) / (12 * dz);

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