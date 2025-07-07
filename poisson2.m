function [p] = poisson2(ni, nj, nk, x, y, z, dt, r, rua, rva, rwa, pa, T, Ta, convu, convv, convw, difu, difv, difw)

    % Inputs:
    %   ni       - Number of grid points in the x-direction  
    %   nj       - Number of grid points in the y-direction  
    %   nk       - Number of grid points in the z-direction  
    %   x        - 1D array of spatial coordinates in the x-direction  
    %   y        - 1D array of spatial coordinates in the y-direction  
    %   z        - 1D array of spatial coordinates in the z-direction  
    %   dt       - Time step size  
    %   r        - 3D array of fluid density:             r    = ρ  
    %   rua      - 3D array of x-direction momentum:      rua  = ρu  
    %   rva      - 3D array of y-direction momentum:      rva  = ρv  
    %   rwa      - 3D array of z-direction momentum:      rwa  = ρw  
    %   pa       - 3D array of pressure from previous time step  
    %   T        - 3D array of temperature at current step  
    %   Ta       - 3D array of temperature at previous step  
    %   convu    - 3D array of x-momentum convective term  
    %   convv    - 3D array of y-momentum convective term  
    %   convw    - 3D array of z-momentum convective term  
    %   difu     - 3D array of x-momentum diffusive term  
    %   difv     - 3D array of y-momentum diffusive term  
    %   difw     - 3D array of z-momentum diffusive term  

    % Author: Jonatan Ismael Eisermann  
    % Date: July 6, 2025.  
    
    p    = pa;

    for i = 2 : ni-1
        for j = 2 : nj-1
            for k = 2 : nk-1
                dx  = 0.5 * (x(i+1) - x(i-1));
                dy  = 0.5 * (y(j+1) - y(j-1));
                dz  = 0.5 * (z(k+1) - z(k-1));
                dx2 = dx * dx;
                dy2 = dy * dy;
                dz2 = dz * dz;

                dconvudx = (convu(i+1,j,k) - convu(i-1,j,k)) / (2 * dx);
                dconvvdy = (convv(i,j+1,k) - convv(i,j-1,k)) / (2 * dy);
                dconvwdz = (convw(i,j,k+1) - convw(i,j,k-1)) / (2 * dz);

                ddifudx  = (difu(i+1,j,k) - difu(i-1,j,k)) / (2 * dx);
                ddifvdy  = (difv(i,j+1,k) - difv(i,j-1,k)) / (2 * dy);
                ddifwdz  = (difw(i,j,k+1) - difw(i,j,k-1)) / (2 * dz);

                px = (pa(i+1,j,k) + pa(i-1,j,k)) / dx2;
                py = (pa(i,j+1,k) + pa(i,j-1,k)) / dy2;
                pz = (pa(i,j,k+1) + pa(i,j,k-1)) / dz2;

                coef     = (dx2 * dy2 * dz2) / (2 * (dx2 * dy2 + dy2 * dz2 + dx2 * dz2));
                
                % continuity equation
                drudx0 = (rua(i+1,j,k) - rua(i-1,j,k)) / (2 * dx);
                drvdy0 = (rva(i,j+1,k) - rva(i,j-1,k)) / (2 * dy);
                drwdz0 = (rwa(i,j,k+1) - rwa(i,j,k-1)) / (2 * dz);

                Da(i,j,k) = drudx0 + drvdy0 + drwdz0; 
                D(i,j,k)  = (r(i,j,k) / T(i,j,k)) * (T(i,j,k) - Ta(i,j,k)) / dt;

                dDdt     = (D(i,j,k) - Da(i,j,k)) / dt;
                convterm = dconvudx + dconvvdy + dconvwdz;
                diffterm = ddifudx + ddifvdy + ddifwdz;
                
                p(i,j,k) = coef * (px + py + pz + dDdt + convterm - diffterm);
            end
        end
    end

end