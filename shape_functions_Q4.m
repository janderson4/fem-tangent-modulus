% ---- Shape function subroutine ----
function [B, detJ] = shape_functions_Q4(xi, eta, coords, nu)
    % Calculates the derivative the the shape functions with respect to the
    % coordinates x and y as well as the integration factor for the 
    % element stiffness.
    % Args:
    %   xi: the xi coordinate of the IP
    %   eta: the eta coordinate of the IP
    %   coords: the coordinates of the 4 nodes of the element, in 
    %       conventional ordering and x,y over the columns (4 x 2)
    % Returns:
    %   B: the derivative of each of the nonzero stress components at this
    %       IP (sigma_x, sigma_y, sigma_xy) with respect to the x and y
    %       coordinates (3 x 8)
    %   detJ: the determinate of the matrix J which is the jacobian of the
    %   x,y coordinates (differentiated with respect to xi, eta)
    
    dNdxi = 0.25 * [  -(1 - eta),  (1 - eta),  (1 + eta), -(1 + eta);
                     -(1 - xi),  -(1 + xi),   (1 + xi),   (1 - xi)];

    % J is dx_dxi (2x2)
    J = dNdxi * coords;
    detJ = det(J);

    dxi_dx=inv(J);

    % since dNdxi is split by row between xi and eta, the transpose of
    % dxi_dx must be used for the chain rule
    dNdx = dxi_dx' * dNdxi;

    % B is the strain-displacement transformation matrix (only including
    % the strains which are calculated from the displacements)
    B = zeros(6,8);
    for k = 1:4
        B(1,2*k-1) = dNdx(1,k);
        B(2,2*k)   = dNdx(2,k);
        B(3,2*k-1) = -nu/(1-nu)*dNdx(2-1,k);
        B(3,2*k)   = -nu/(1-nu)*dNdx(2,k);
        B(6,2*k-1) = dNdx(2,k);
        B(6,2*k)   = dNdx(1,k);
    end
end