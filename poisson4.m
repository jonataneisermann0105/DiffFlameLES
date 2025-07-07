function [p] = poisson4(ni, nj, nk, x, y, z, dt, r, rua, rva, rwa, pa, T, Ta, convu, convv, convw, difu, difv, difw)

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

    % Outputs:
    %   p        - 3D array of updated pressure field  

    % Author: Jonatan Ismael Eisermann  
    % Date: July 6, 2025.  

    p    = pa;

    for i = 2 : ni-1
        for j = 2 : nj-1
            for k = 2 : nk-1

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

                dx  = (x(i+1) - x(i-1)) / 2;
                dy  = (y(j+1) - y(j-1)) / 2;
                dz  = (z(k+1) - z(k-1)) / 2;
                dx2 = dx * dx;
                dy2 = dy * dy;
                dz2 = dz * dz;

                dconvudx = (- convu(ie,j,k) + 8 * convu(i+1,j,k) - 8 * convu(i-1,j,k) + convu(iw,j,k)) / (12 * dx);
                dconvvdy = (- convv(i,jn,k) + 8 * convv(i,j+1,k) - 8 * convv(i,j-1,k) + convv(i,js,k)) / (12 * dy);
                dconvwdz = (- convw(i,j,kf) + 8 * convw(i,j,k+1) - 8 * convw(i,j,k-1) + convw(i,j,kb)) / (12 * dz);

                ddifudx = (- difu(ie,j,k) + 8 * difu(i+1,j,k) - 8 * difu(i-1,j,k) + difu(iw,j,k)) / (12 * dx);
                ddifvdy = (- difv(i,jn,k) + 8 * difv(i,j+1,k) - 8 * difv(i,j-1,k) + difv(i,js,k)) / (12 * dy);
                ddifwdz = (- difw(i,j,kf) + 8 * difw(i,j,k+1) - 8 * difw(i,j,k-1) + difw(i,j,kb)) / (12 * dz);

                px   = (- pa(ie,j,k) + 16 * pa(i+1,j,k) + 16 * pa(i-1,j,k) - pa(iw,j,k)) / (12 * dx2);
                py   = (- pa(i,jn,k) + 16 * pa(i,j+1,k) + 16 * pa(i,j-1,k) - pa(i,js,k)) / (12 * dy2);
                pz   = (- pa(i,j,kf) + 16 * pa(i,j,k+1) + 16 * pa(i,j,k-1) - pa(i,j,kb)) / (12 * dz2);

                coef     = (dx2 * dy2 * dz2) / (5/2 * (dx2 * dy2 + dy2 * dz2 + dx2 * dz2));
                
                 % continuity equation
                drudx0 = (- rua(ie,j,k) + 8 * rua(i+1,j,k) - 8 * rua(i-1,j,k) + rua(iw,j,k)) / (12 * dx);
                drvdy0 = (- rva(i,jn,k) + 8 * rva(i,j+1,k) - 8 * rva(i,j-1,k) + rva(i,js,k)) / (12 * dy);
                drwdz0 = (- rwa(i,j,kf) + 8 * rwa(i,j,k+1) - 8 * rwa(i,j,k-1) + rwa(i,j,kb)) / (12 * dz);

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