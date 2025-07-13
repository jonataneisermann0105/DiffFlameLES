function [difu,difv,difw,difZ] = diffusive2(ni, nj, nk, x, y, z, ra, rua, rva, rwa, rZa, T, Tin, Cs, Sc, rout, nuout)

    % Inputs:
    %   ni    - Number of grid points in the x-direction
    %   nj    - Number of grid points in the y-direction
    %   nk    - Number of grid points in the z-direction
    %   x     - 1D array of spatial coordinates in the x-direction
    %   y     - 1D array of spatial coordinates in the y-direction
    %   z     - 1D array of spatial coordinates in the z-direction
    %   ra    - 3D array of fluid density:         ra  = ρ
    %   rua   - 3D array of momentum in x:         rua = ρu
    %   rva   - 3D array of momentum in y:         rva = ρv
    %   rwa   - 3D array of momentum in z:         rwa = ρw
    %   rZa   - 3D array of a conserved scalar:    rZa = ρZ
    %   T     - 3D array of temperature field
    %   Tin   - Initial temperature
    %   Cs    - Smagorinsky constant
    %   Sc    - Schmidt number
    %   rout  - Reference density
    %   nuout - Reference kinematic viscosity

    % Outputs:
    %   difu  - Diffusive term of the x-momentum equation
    %   difv  - Diffusive term of the y-momentum equation
    %   difw  - Diffusive term of the z-momentum equation
    %   difZ  - Diffusive term of the mixture fraction equation

    % Author: Jonatan Ismael Eisermann
    % Date: July 6, 2025.

    tau_xx = zeros(ni,nj,nk);
    tau_yy = zeros(ni,nj,nk);
    tau_zz = zeros(ni,nj,nk);
    tau_xy = zeros(ni,nj,nk);
    tau_xz = zeros(ni,nj,nk);
    tau_yz = zeros(ni,nj,nk);
    tau_Zx = zeros(ni,nj,nk);
    tau_Zy = zeros(ni,nj,nk);
    tau_Zz = zeros(ni,nj,nk);

    for k = 2 : nk-1
        for j = 2 : nj-1
            for i = 2 : ni-1
                
                dx  = (x(i+1) - x(i-1)) / 2;
                dy  = (y(j+1) - y(j-1)) / 2;
                dz  = (z(k+1) - z(k-1)) / 2;
                Delta = (dx * dy * dz) ^ (1/3);

                dudx = (rua(i+1,j,k) / ra(i+1,j,k) - rua(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);
                dvdx = (rva(i+1,j,k) / ra(i+1,j,k) - rva(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);
                dwdx = (rwa(i+1,j,k) / ra(i+1,j,k) - rwa(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);
                dZdx = (rZa(i+1,j,k) / ra(i+1,j,k) - rZa(i-1,j,k) / ra(i-1,j,k)) / (2 * dx);

                dudy = (rua(i,j+1,k) / ra(i,j+1,k) - rua(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);
                dvdy = (rva(i,j+1,k) / ra(i,j+1,k) - rva(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);
                dwdy = (rwa(i,j+1,k) / ra(i,j+1,k) - rwa(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);
                dZdy = (rZa(i,j+1,k) / ra(i,j+1,k) - rZa(i,j-1,k) / ra(i,j-1,k)) / (2 * dy);

                dudz = (rua(i,j,k+1) / ra(i,j,k+1) - rua(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);  
                dvdz = (rva(i,j,k+1) / ra(i,j,k+1) - rva(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);                        
                dwdz = (rwa(i,j,k+1) / ra(i,j,k+1) - rwa(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);  
                dZdz = (rZa(i,j,k+1) / ra(i,j,k+1) - rZa(i,j,k-1) / ra(i,j,k-1)) / (2 * dz);

                nut   = (Cs * Delta)^2 * sqrt(2*(dudx^2 + (dudy^2)/2 + dudy * dvdx + (dvdx^2)/2 + dvdy^2 + (dudz^2)/2 + dudz * dwdx + (dwdx^2)/2 + (dvdz^2)/2 + dvdz * dwdy + (dwdy^2)/2 + dwdz^2));
                nul   = ((rout * nuout) / ra(i,j,k)) * (T(i,j,k) / Tin)^(0.7);
                nu    = nul + nut;

                tau_xx(i,j,k) = ra(i,j,k) * nu * (2 * dudx - (2/3) * (dudx + dvdy + dwdz));
                tau_yy(i,j,k) = ra(i,j,k) * nu * (2 * dvdy - (2/3) * (dudx + dvdy + dwdz));
                tau_zz(i,j,k) = ra(i,j,k) * nu * (2 * dwdz - (2/3) * (dudx + dvdy + dwdz));

                tau_xy(i,j,k) = ra(i,j,k) * nu * (dudy + dvdx);
                tau_xz(i,j,k) = ra(i,j,k) * nu * (dudz + dwdx); 
                tau_yz(i,j,k) = ra(i,j,k) * nu * (dvdz + dwdy);

                tau_Zx(i,j,k)  = ra(i,j,k) * nu / Sc * dZdx;
                tau_Zy(i,j,k)  = ra(i,j,k) * nu / Sc * dZdy;
                tau_Zz(i,j,k)  = ra(i,j,k) * nu / Sc * dZdz;

            end
        end
    end

    for k = 1 : nk
        for j = 1 : nj
            tau_xx(ni,j,k) = tau_xx(ni-1,j,k);
            tau_yy(ni,j,k) = tau_yy(ni-1,j,k);
            tau_zz(ni,j,k) = tau_zz(ni-1,j,k);
            tau_xy(ni,j,k) = tau_xy(ni-1,j,k);
            tau_xz(ni,j,k) = tau_xz(ni-1,j,k);
            tau_yz(ni,j,k) = tau_yz(ni-1,j,k);
            tau_Zx(ni,j,k) = tau_Zx(ni-1,j,k);
            tau_Zy(ni,j,k) = tau_Zy(ni-1,j,k);
            tau_Zz(ni,j,k) = tau_Zz(ni-1,j,k);
        end
    end

    difu = zeros(ni,nj,nk);
    difv = zeros(ni,nj,nk);
    difw = zeros(ni,nj,nk);
    difZ = zeros(ni,nj,nk);

    for k = 2 : nk-1
        for j = 2 : nj-1
            for i = 2 : ni-1

                dtau_xxdx = (tau_xx(i+1,j,k) - tau_xx(i-1,j,k)) / (2 * dx);
                dtau_xydx = (tau_xy(i+1,j,k) - tau_xy(i-1,j,k)) / (2 * dx);
                dtau_xzdx = (tau_xz(i+1,j,k) - tau_xz(i-1,j,k)) / (2 * dx);
                dtau_Zxdx = (tau_Zx(i+1,j,k) - tau_Zx(i-1,j,k)) / (2 * dx);
                    
                dtau_xydy = (tau_xy(i,j+1,k) - tau_xy(i,j-1,k)) / (2 * dy);
                dtau_yydy = (tau_yy(i,j+1,k) - tau_yy(i,j-1,k)) / (2 * dy);
                dtau_yzdy = (tau_yz(i,j+1,k) - tau_yz(i,j-1,k)) / (2 * dy);
                dtau_Zydy = (tau_Zy(i,j+1,k) - tau_Zy(i,j-1,k)) / (2 * dy);

                dtau_xzdz = (tau_xz(i,j,k+1) - tau_xz(i,j,k-1)) / (2 * dz);
                dtau_yzdz = (tau_yz(i,j,k+1) - tau_yz(i,j,k-1)) / (2 * dz);
                dtau_zzdz = (tau_zz(i,j,k+1) - tau_zz(i,j,k-1)) / (2 * dz);
                dtau_Zzdz = (tau_Zz(i,j,k+1) - tau_Zz(i,j,k-1)) / (2 * dz);
                    
                difu(i,j,k) = dtau_xxdx + dtau_xydy + dtau_xzdz;
                difv(i,j,k) = dtau_xydx + dtau_yydy + dtau_yzdz; 
                difw(i,j,k) = dtau_xzdx + dtau_yzdy + dtau_zzdz; 
                difZ(i,j,k) = dtau_Zxdx + dtau_Zydy + dtau_Zzdz;

            end
        end
    end

    for k = 1 : nk
        for j = 1 : nj
            difu(ni,j,k) = difu(ni-1,j,k);
            difv(ni,j,k) = difv(ni-1,j,k);
            difw(ni,j,k) = difw(ni-1,j,k);
            difZ(ni,j,k) = difZ(ni-1,j,k);
        end
    end

end