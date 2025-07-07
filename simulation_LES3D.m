%*************************************************************************%
       %% COMPUTATIONAL SIMULATION OF THE SANDIA DME D FLAME %%
%*************************************************************************%

% Description:
% This MATLAB code simulates a turbulent diffusion flame at low Mach 
% number using Large Eddy Simulation (LES). It employs a finite difference 
% method for spatial discretization and a projection-type algorithm for 
% temporal integration. Subgrid-scale turbulence effects are modeled using 
% the classical Smagorinsky model.

% Notes:
% Make sure all auxiliary functions listed below are available in the 
% working directory, as each one is responsible for a specific numerical 
% computation or post-processing task:
% - convective2.m: computes convective terms using second-order central 
%                  finite differences;
% - convective4.m: computes convective terms using fourth-order central 
%                  finite differences;
% - pressure2.m:   calculates the pressure gradient using second-order 
%                  finite differences;
% - pressure4.m:   calculates the pressure gradient using fourth-order 
%                  finite differences;
% - diffusive2.m:  computes viscous/diffusive terms using second-order 
%                  finite differences;
% - diffusive4.m:  computes viscous/diffusive terms using fourth-order 
%                  finite differences;
% - poisson2.m:    solves the Poisson equation for pressure correction 
%                  using second-order discretization;
% - poisson4.m:    solves the Poisson equation for pressure correction 
%                  using fourth-order discretization;
% - tempdens.m:    computes temperature and density fields based on the 
%                  mixture fraction and thermodynamic principles;
% - animation.m:   generates an animated .gif showing the temporal 
%                  evolution of key flow variables (velocity magnitude, 
%                  pressure, mixture fraction, temperature).

% Author: Jonatan Ismael Eisermann  
% Date: July 6, 2025. 

clear; clc;                  % data cleaning

disp('         COMPUTER SIMULATION OF THE SANDIA DME D FLAME           ');
disp(' ');
disp('Computer code developed by Jonatan Ismael Eisermann as part of his');
disp('doctoral research at the Federal University of Rio Grande do Sul, under');
disp('the supervision of Prof. Álvaro Luiz de Bortoli, which includes a');
disp('research period at the University of Seville under the supervision of');
disp('Prof. Samuele Rubino.');
disp(' ');
disp('This code uses Large Eddy Simulation (LES) with the Smagorinsky model to');
disp('solve the compressible conservation equations for a turbulent reactive');
disp('flow. The equations are formulated in conservative and dimensionless form,');
disp('assuming a low Mach number regime, where pressure variations are minimal.');
disp('This assumption enables the use of simplified numerical methods, such as');
disp('the projective-type methods.');
disp(' ');
disp('Computater code developed in MATLAB, version R2024b.');
disp(' ');

s_order = input('Please enter the desired spatial discretization order (2 or 4): ');
disp(' ');

while ~(s_order == 2 || s_order == 4)
    disp('Invalid input. Please enter 2 or 4.');
    s_order = input('Please enter the desired spatial discretization order (2 or 4): ');
end

fprintf('You have selected a spatial discretization order of: %d\n', s_order);

t_order = input('Please enter the desired temporal integration order (1 or 2): ');
disp(' ');

while ~(t_order == 1 || t_order == 2)
    disp('Invalid input. Please enter 1 or 2.');
    t_order = input('Please enter the desired temporal integration order (1 or 2): ');
end

fprintf('You have selected a temporal integration order of: %d\n', t_order);

itmax = input('Please enter the maximum number of iterations: ');
disp(' ');

while ~(isnumeric(itmax) && itmax > 0 && floor(itmax) == itmax)
    disp('Invalid input. Please enter a positive integer.');
    itmax = input('Please enter the maximum number of iterations: ');
end

ts = input('Please enter the dimensional time step size: ');
disp(' ');

while ~(isnumeric(ts) && ts > 0)
    disp('Invalid input. Please enter a positive number.');
    ts = input('Please enter the dimensional time step size: ');
end

%*************************************************************************%
% Characteristic scales
%*************************************************************************%
lcar  = 0.00745;                   % main nozzle diameter (m)
ucar  = 45.9;                      % fuel injection velocity (m/s)
rcar  = 1.323;                     % fuel density (kg/m^3)
pcar  = rcar * ucar^2;             % pressure (Pa)
R     = 101325 / (1.184 * 293.15); % R = p / (r * T);
Tcar  = ucar^2 / R;                % temperature (K)
nu0   = 0.00001167;                % fuel kinematic viscosity (m^2/s)

%*************************************************************************%
% Dimensionless quantities
%*************************************************************************%
Sc    = 0.7;                 % Schmidt number
uin   = 0.9 / ucar;          % velocity u at the air inlet
vin   = 0 / ucar;			 % inlet velocity v  
win   = 0 / ucar;            % inlet velocity w
uout  = 45.9 / ucar;         % velocity u at the fuel inlet
upil  = 1.1 / ucar;          % velocity u at the pilot inlet
rin   = 1.184 / rcar;        % density at the air and pilot inlet
rout  = 1.323 / rcar;        % density at the fuel inlet 
nuout = nu0 / (ucar * lcar); % kinematic viscosity
pin   = 101325 / pcar;       % inlet pressure 
Tin   = 293.15 / Tcar;       % inlet temperature
dt    = ts * ucar / lcar;    % time step
Zst   = 0.35;                % stoichiometric mixture fraction
Tad   = 2100 / Tcar;         % adiabatic flame temperature
heat  = Tad - Tin;           % Maximum temperature increase due to reaction  

%*************************************************************************%
% Spatial domain
%*************************************************************************%
comp  = 50;     	         % length, in diameters d, in the x direction
alt   = 8;	 		         % length, in diameters d, in the y direction  
esp   = alt;                 % length, in diameters d, in the z direction
ri    = (lcar / lcar) / 2;   % fuel nozzle radius
rf    = (0.0182 / lcar) / 2; % pilot annulus radius
xf    = comp / 50;           % length of the fuel injection cylinder

%*************************************************************************%
% Computational mesh
%*************************************************************************%
ni    = 99;                  % number of mesh points in the x direction
nj    = 99;                  % number of mesh points in the y direction
nk    = nj;                  % number of mesh points in the z direction
fati  = 1;                   % mesh concentration factor in the x direction
fatj  = 1;                   % mesh concentration factor in the y direction
fatk  = 1;                   % mesh concentration factor in the z direction

veti(1) = 1;
somai   = 0;
for i = 2:ni
    veti(i) = veti(i-1) * fati;
    if(i < ni/2)
        veti(i) = veti(i-1);
    end
    somai = somai + veti(i);
end
x(1) = 0;
for i = 2:ni
    x(i) = x(i-1) + comp * veti(i) / somai;
end

vetj(1) = 1;
somaj   = 0;
for j = 2:nj
    vetj(j) = vetj(j-1) * fatj;
    if(j > nj/2)
        vetj(j) = vetj(nj-j+1);
    end
    somaj = somaj + vetj(j);
end
y(1) = 0;
for j = 2:nj
    y(j) = y(j-1) + alt * vetj(j) / somaj;
end

vetk(1) = 1;
somak   = 0;
for k = 2:nk
    vetk(k) = vetk(k-1) * fatk;
    if(k > nk/2)
        vetk(k) = vetk(nk-k+1);
    end
    somak = somak + vetk(k);
end
z(1) = 0;
for k = 2:nk
    z(k) = z(k-1) + esp * vetk(k) / somak;
end

%*************************************************************************%
% Initial conditions
%*************************************************************************%
for k = 1 : nk
    for j = 1 : nj
        for i = 1 : ni
            raio           = sqrt((y(j) - alt / 2)^2 + (z(k) - esp / 2)^2);
            pa(i,j,k)      = pin;
            ra(i,j,k)      = rin;
            rua(i,j,k)     = 0;
            rva(i,j,k)     = rin * vin;
            rwa(i,j,k)     = rin * win;
            rZa(i,j,k)     = 0;
            Ta(i,j,k)      = Tin;

            if x(i) <= xf
                rua(i,j,k) = rin * uin;
                if raio < rf
                    rua(i,j,k) = rin * upil;
                    rZa(i,j,k) = rin * 0.27;
                    if raio < ri
                        ra(i,j,k)  = rout;
                        rua(i,j,k) = rout * uout * (1 - (raio / ri)^2);
                        rZa(i,j,k) = rout * 1; 
                    end
                end
            end

            p(i,j,k)      = pa(i,j,k);
            r(i,j,k)      = ra(i,j,k);
            ru(i,j,k)     = rua(i,j,k);
            rv(i,j,k)     = rva(i,j,k);
            rw(i,j,k)     = rwa(i,j,k);
            rZ(i,j,k)     = rZa(i,j,k);
            T(i,j,k)      = Ta(i,j,k);
        end
    end
end

%*************************************************************************%
% Temporal numerical integration
%*************************************************************************%
Cs      = 0.3;                  % Smagorisnky coefficient
beta    = 0.999;                % relaxation factor for density
h       = figure(1);            % figure for the GIF animation

for it = 1 : itmax
    time = dt * it;
    res  = 0;

    if s_order == 2
    
        % Convective terms
        [convu, convv, convw, convZ] = convective2(ni, nj, nk, x, y, z, ra, rua, rva, rwa, rZa);
    
        % Pressure gradient
        [dpdx, dpdy, dpdz] = pressure2(ni, nj, nk, x, y, z, pa);

        % Diffusive terms
        [difu, difv, difw, difZ] = diffusive2(ni, nj, nk, x, y, z, ra, rua, rva, rwa, rZa, T, Tin, Cs, Sc, rout, nuout);

    elseif s_order == 4

        % Convective terms
        [convu, convv, convw, convZ] = convective4(ni, nj, nk, x, y, z, ra, rua, rva, rwa, rZa);
    
        % Pressure gradient
        [dpdx, dpdy, dpdz] = pressure4(ni, nj, nk, x, y, z, pa);

        % Diffusive terms
        [difu, difv, difw, difZ] = diffusive4(ni, nj, nk, x, y, z, ra, rua, rva, rwa, rZa, Ta, Tin, Cs, Sc, rout, nuout);

    end 

    for k = 2 : nk-1
        for j = 2 : nj-1
            for i = 2 : ni-1

                fu(i,j,k) = - dpdx(i,j,k) + difu(i,j,k) - convu(i,j,k);
                fv(i,j,k) = - dpdy(i,j,k) + difv(i,j,k) - convv(i,j,k);
                fw(i,j,k) = - dpdz(i,j,k) + difw(i,j,k) - convw(i,j,k);
                fZ(i,j,k) =                 difZ(i,j,k) - convZ(i,j,k);

                ru(i,j,k) = rua(i,j,k) + dt * fu(i,j,k);
                rv(i,j,k) = rva(i,j,k) + dt * fv(i,j,k);
                rw(i,j,k) = rwa(i,j,k) + dt * fw(i,j,k);
                rZ(i,j,k) = rZa(i,j,k) + dt * fZ(i,j,k);

                if rZ(i,j,k) < 0
                    rZ(i,j,k) = 0;
                end
    
                if (rZ(i,j,k) / ra(i,j,k)) > 1
                    rZ(i,j,k) = 1 * ra(i,j,k);
                end

            end
        end
    end

    % Calculation of temperature and density
    [Z, T, r] = tempdens(alt, esp, ni, nj, nk, x, y, z, ra, rZ, Ta, Zst, Tin, heat, ri, xf, beta, pin, rout);
    
    if t_order == 2

        if s_order == 2
    
            % Convective terms
            [convu2, convv2, convw2, convZ2] = convective2(ni, nj, nk, x, y, z, r, ru, rv, rw, rZ);

            % Diffusive terms
            [difu2, difv2, difw2, difZ2] = diffusive2(ni, nj, nk, x, y, z, r, ru, rv, rw, rZ, T, Tin, Cs, Sc, rout, nuout);

        elseif s_order == 4

            % Convective terms
            [convu2, convv2, convw2, convZ2] = convective4(ni, nj, nk, x, y, z, r, ru, rv, rw, rZ);

            % Diffusive terms
            [difu2, difv2, difw2, difZ2] = diffusive4(ni, nj, nk, x, y, z, r, ru, rv, rw, rZ, T, Tin, Cs, Sc, rout, nuout);

        end 

        convu = (convu + convu2) / 2;
        convv = (convv + convv2) / 2;
        convw = (convw + convw2) / 2;
        convZ = (convZ + convZ2) / 2;
        
        difu  = ( difu + difu2 ) / 2;
        difv  = ( difv + difv2 ) / 2;
        difw  = ( difw + difw2 ) / 2;
        difZ  = ( difZ + difZ2 ) / 2;

        for k = 2 : nk-1
            for j = 2 : nj-1
                for i = 2 : ni-1

                    fu2(i,j,k) = - dpdx(i,j,k) + difu2(i,j,k) - convu2(i,j,k);
                    fv2(i,j,k) = - dpdy(i,j,k) + difv2(i,j,k) - convv2(i,j,k);
                    fw2(i,j,k) = - dpdz(i,j,k) + difw2(i,j,k) - convw2(i,j,k);
                    fZ2(i,j,k) =                 difZ2(i,j,k) - convZ2(i,j,k);

                    ru(i,j,k)  = rua(i,j,k) + 0.5 * dt * (fu(i,j,k) + fu2(i,j,k));
                    rv(i,j,k)  = rva(i,j,k) + 0.5 * dt * (fv(i,j,k) + fv2(i,j,k));
                    rw(i,j,k)  = rwa(i,j,k) + 0.5 * dt * (fw(i,j,k) + fw2(i,j,k));
                    rZ(i,j,k)  = rZa(i,j,k) + 0.5 * dt * (fZ(i,j,k) + fZ2(i,j,k));

                end
            end
        end

    end

    % Low Mach number continuity residual
    for k = 2 : nk-1
        for j = 2 : nj-1
            for i = 2 : ni-1

                dx  = 0.5 * (x(i+1) - x(i-1));
                dy  = 0.5 * (y(j+1) - y(j-1));
                dz  = 0.5 * (z(k+1) - z(k-1));

                drudx = (ru(i+1,j,k) - ru(i-1,j,k)) / (2 * dx);
                drvdy = (rv(i,j+1,k) - rv(i,j-1,k)) / (2 * dy);
                drwdz = (rw(i,j,k+1) - rw(i,j,k-1)) / (2 * dz);
                dTdt  = (T(i,j,k) - Ta(i,j,k)) / dt;

                rest = drudx + drvdy + drwdz - (r(i,j,k) / T(i,j,k)) * dTdt; 

                if rest > res 
                    res = rest;
                end

            end
        end
    end

    % pressure calculation using a Poisson equation
    if s_order == 2
        [p] = poisson2(ni, nj, nk, x, y, z, dt, r, rua, rva, rwa, pa, T, Ta, convu, convv, convw, difu, difv, difw);
    elseif s_order ==4 
        [p] = poisson4(ni, nj, nk, x, y, z, dt, r, rua, rva, rwa, pa, T, Ta, convu, convv, convw, difu, difv, difw);
    end
    
    %*********************************************************************%
    % Boundary conditions
    %*********************************************************************%
    % (0,j,k,t) to (comp,j,k,t)
    for i = 1 : ni
        for j = 1 : nj
            for k = 1 : nk
                raio = sqrt((y(j)-alt/2)^2 + (z(k)-esp/2)^2);
                if x(i) <= xf
                    p(1,j,k)      = p(2,j,k);
                    r(i,j,k)      = rin;
                    ru(i,j,k)     = rin * uin;
                    rv(i,j,k)     = rin * vin;
                    rw(i,j,k)     = rin * win;
                    rZ(i,j,k)     = rin * 0;

                    if raio < rf
                        ru(i,j,k) = rin * upil;
                        rZ(i,j,k) = rin * 0.27;
                        if raio < ri
                            r(i,j,k)  = rout;
                            ru(i,j,k) = rout * uout * (1 - (raio/ri)^2);
                            rZ(i,j,k) = rout * 1;
                        end
                    end
                end
            end
        end
    end

    % (comp,y,z,t)
    for j = 1 : nj
        for k = 1 : nk
            raio           = sqrt((y(j)-alt/2)^2 + (z(k)-esp/2)^2);
            p(ni,j,k)      = pin;
            r(ni,j,k)      = r(ni-1,j,k);
            ru(ni,j,k)     = ru(ni-1,j,k);                       
            rv(ni,j,k)     = rv(ni-1,j,k);
            rw(ni,j,k)     = rw(ni-1,j,k);
            rZ(ni,j,k)     = rZ(ni-1,j,k);
        end
    end

    % (x,0,z,t) and (x,alt,z,t)
    for k = 1 : nk
        for i = 1 : ni
            p(i,1,k)      = p(i,2,k);
            p(i,nj,k)     = p(i,nj-1,k);
            r(i,1,k)      = r(i,2,k);
            r(i,nj,k)     = r(i,nj-1,k);
            ru(i,1,k)     = 0; %ru(i,2,k);
            ru(i,nj,k)    = 0; ru(i,nj-1,k);
            rv(i,1,k)     = 0; %rv(i,2,k);
            rv(i,nj,k)    = 0; %rv(i,nj-1,k);
            rw(i,1,k)     = 0; %rw(i,2,k);
            rw(i,nj,k)    = 0; %rw(i,nj-1,k);
            rZ(i,1,k)     = rZ(i,2,k);
            rZ(i,nj,k)    = rZ(i,nj-1,k);
        end
    end

    % (x,y,0,t) and (x,y,esp,t)
    for j = 1 : nj
        for i = 1 : ni
            p(i,j,1)      = p(i,j,2);
            p(i,j,nk)     = p(i,j,nk-1);
            r(i,j,1)      = r(i,j,2);
            r(i,j,nk)     = r(i,j,nk-1);
            ru(i,j,1)     = 0; %ru(i,j,2);
            ru(i,j,nk)    = 0; %ru(i,j,nk-1);
            rv(i,j,1)     = 0; %rv(i,j,2);
            rv(i,j,nk)    = 0; %rv(i,j,nk-1);
            rw(i,j,1)     = 0; %rw(i,j,2);
            rw(i,j,nk)    = 0; %rw(i,j,nk-1);
            rZ(i,j,1)     = rZ(i,j,2);
            rZ(i,j,nk)    = rZ(i,j,nk-1);
        end
    end

    for k = 1 : nk
        for j = 1 : nj
            for i = 1 : ni
                ra(i,j,k)      = r(i,j,k);
                pa(i,j,k)      = p(i,j,k);
                rua(i,j,k)     = ru(i,j,k);
                rva(i,j,k)     = rv(i,j,k);
                rwa(i,j,k)     = rw(i,j,k);
                rZa(i,j,k)     = rZ(i,j,k);
                Ta(i,j,k)      = T(i,j,k);
            end
        end
    end

    % Screen output: iteration number and low Mach continuity residual
    if mod(it, 1) == 0
        disp(['it ', num2str(it), '     Continuity residual ', num2str(res)]);
    end

    kmeio = (nk+1) / 2;
    for i = 1 : ni
        for j = 1 : nj
            for k = 1 : nk
                u(i, j, k)   = ru(i, j, k) / r(i, j, k);
                v(i, j, k)   = rv(i, j, k) / r(i, j, k);
                w(i, j, k)   = rw(i, j, k) / r(i, j, k);
                vel(i, j, k) = sqrt((u(i, j, k))^2 + (v(i, j, k))^2 + (w(i, j, k))^2);
            end
            vel2d(i,j) = sqrt((u(i, j, kmeio))^2 + (v(i, j, kmeio))^2 + (w(i, j, kmeio))^2);
            p2d(i,j)    = p(i, j , kmeio);
            Z2d(i,j)    = Z(i, j, kmeio);
            T2d(i,j)    = T(i, j, kmeio);
            r2d(i,j)    = r(i, j, kmeio);
        end
    end
 
    % Generate the plots for the animation

    if (mod(it, 100) == 0) || (it == 1)
        % Plot the velocity magnitude
        cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
        subplot(4, 1, 1);
        contourf(x, y, vel2d', 23, 'LineColor', 'none');
        titulo0 = sprintf('Dimensionless velocity at t*=%.4f', time);
        title(titulo0);
        axis('equal',[0 comp 0 alt]);
        xlabel('x/d'); ylabel('y/d');
        colorbar;

        % Plot the pressure
        cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
        subplot(4, 1, 2);
        contourf(x, y, p2d', 23, 'LineColor', 'none');
        titulo = sprintf('Dimensionless pressure at t*=%.4f', time);
        title(titulo);
        axis('equal', [0 comp 0 alt]);
        xlabel('x/d'); ylabel('y/d');
        colorbar;

         % Plot the mixture fraction
        cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
        subplot(4, 1, 3);
        contourf(x, y, Z2d', 23, 'LineColor', 'none');
        titulo0 = sprintf('Mixture fraction at t*=%.4f', time);
        title(titulo0);
        axis('equal',[0 comp 0 alt]);
        xlabel('x/d'); ylabel('y/d');
        colorbar;

         % Plot the temperature
        cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
        subplot(4, 1, 4);
        contourf(x, y, T2d', 23, 'LineColor', 'none');
        titulo0 = sprintf('Dimensionless temperature at t*=%.4f', time);
        title(titulo0);
        axis('equal',[0 comp 0 alt]);
        xlabel('x/d'); ylabel('y/d');
        colorbar;
        drawnow; % Atualiza a figura
        animation(h, it, 'Simulation.gif');
        pause(0.000005);
    end

end		