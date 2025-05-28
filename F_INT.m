function [f_int, A, SIG, ALPHA, R]=F_INT(SIG,u_diff,Ce,ALPHA,R,mu,H,beta,nodes,elements,t,nu)
    % Computes the internal forces in the structure and the algorithmic
    % consistent tangent operator for given internal state variables and 
    % change in displacement. 
    % Args:
    %   SIG: the stress at all IPs before this step (6 x nIP)
    %   u_diff: the trial change in displacement at all dofs during this 
    %       step (2*ndof x 1)
    %   Ce: elastic tangent operator
    %   ALPHA: back stress for all IPs before this step (6 x nIP)
    %   R: radius for all IPs before this step (1 x nIP)
    %   mu: shear modulus
    %   H: hardening modulus
    %   beta: kinematic-isotropic interpolation factor (1=isotropic)
    %   nodes: the nodal coordinates (nnode x 2)
    %   elements: the node numbers for each element (nel x 4)
    %   t: the plate thickness
    % Returns:
    %   f_int: the vector of internal forces (ndof x 1)
    %   A: the algorithmic consistent tangent operator (ndof x ndof)
    %   SIG: the stress at all IPs after this step (6 x nIP)
    %   ALPHA: back stress for all IPs after this step (6 x nIP)
    %   R: radius for all IPs after this step (1 x nIP)

    nel=size(elements,1);
    ndof=size(u_diff,1);

    A = sparse(ndof, ndof);
    f_int = zeros(ndof,1);

    SIG0=SIG;
    ALPHA0=ALPHA;
    R0=R;
    
    for e = 1:nel
        Ae = zeros(8,8);
        f_inte=zeros(8,1);
        coords = nodes(elements(e,:), :);

        % get the dof numbers for this element in node order then x,y order
        edof = reshape([2*elements(e,:)-1; 2*elements(e,:)], 1, []);

        % the change in displacements at this element (8x1)
        u_el=u_diff(edof,:);
        
        for IP = 1:4
            IP_global=(e-1)*4+IP; % the global integration point number
            [xi, eta]=gauss(IP);
            [B, detJ] = shape_functions_Q4(xi, eta, coords, nu);

            % the complete strain vector (6 x 1)
            EPSI_el=B*u_el;

            [SIG(:,IP_global), ALPHA(:,IP_global), R(:,IP_global), C_alg] =...
                sigma(SIG0(:,IP_global),Ce,EPSI_el,ALPHA0(:,IP_global),R0(:,IP_global),mu,H,beta);

            Ae = Ae + B' * C_alg * B * detJ * t;
            f_inte = f_inte + B' * SIG(:,IP_global) * detJ * t;
        end
        
        % Assembly into global
        A(edof, edof) = A(edof, edof) + Ae;
        f_int(edof)=f_int(edof)+f_inte;
    end
end