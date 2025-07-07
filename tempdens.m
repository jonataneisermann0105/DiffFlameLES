function [Z, T, r] = tempdens(alt, esp, ni, nj, nk, x, y, z, ra, rZ, Ta, Zst, Tin, heat, ri, xf, beta, pin, rout)

    % Inputs:
    %   alt     - Domain height in the y-direction  
    %   esp     - Domain width in the z-direction 
    %   ni      - Number of grid points in the x-direction  
    %   nj      - Number of grid points in the y-direction  
    %   nk      - Number of grid points in the z-direction  
    %   x       - 1D array of spatial coordinates in the x-direction  
    %   y       - 1D array of spatial coordinates in the y-direction  
    %   z       - 1D array of spatial coordinates in the z-direction  
    %   ra      - 3D array of density from the previous time step (ρ)  
    %   rZ      - 3D array of mixture fraction from the previous step (ρZ)  
    %   Ta      - 3D array of temperature from the previous time step  
    %   Zst     - Stoichiometric mixture fraction  
    %   Tin     - Inlet temperature  
    %   heat    - Maximum temperature increase due to reaction  
    %   ri      - Burner inner radius (radial cutoff for flame core)  
    %   xf      - Axial position of flame front (used for core cutoff)  
    %   beta    - Relaxation factor (0 ≤ beta ≤ 1) for temperature update  
    %   pin     - Inlet pressure  
    %   rout    - Assigned density value inside the burner core  
 
    % Outputs:
    %   T       - 3D array of updated temperature field  
    %   r       - 3D array of updated density field  
    %   Z       - 3D array of updated mixture fraction field  

    % Author: Jonatan Ismael Eisermann  
    % Date: July 6, 2025.  

    for k = 1 : nk
        for j = 1 : nj
            for i = 1 : ni
                raio     = sqrt((y(j)-alt/2)^2 + (z(k)-esp/2)^2);
                Z(i,j,k) = rZ(i,j,k) / ra(i,j,k);

	            if Z(i,j,k) > 1 
                    Z(i,j,k) = 1;
                end    

                if Z(i,j,k) >= Zst
                    T(i,j,k) = Tin + heat*(((1 - Z(i,j,k)) / (1 - Zst)));
                end

                if Z(i,j,k) < Zst
                    T(i,j,k) = Tin + heat*((Z(i,j,k) / Zst));
                end

	            if (raio < ri) && (x(i)<xf)
                    T(i,j,k) = Tin;
                end 

                % Relaxation for temperature
                T(i,j,k)  = beta * Ta(i,j,k) + (1 - beta) * T(i,j,k);

                % Calculation of the specific mass based on the steady equation
                r(i,j,k) = pin / T(i,j,k); 

                if x(i) <= xf
                        
                    if raio < ri 
                        r(i,j,k) = rout;
                    end

                end
                
            end
        end
    end

end