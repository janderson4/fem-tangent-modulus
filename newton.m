function [u, SIG, ALPHA, R] = newton(SIG,u,Ce,ALPHA,R,mu,H,beta,f_ext,nodes,elements,t,nu,free_dofs)
    % Uses Newton-Raphson iterations to compute the displacement and new state
    % variables after a loading or unloading step, using a combined
    % kinematic-isotropic hardening model
    % Args:
    %   SIG: the stress at all IPs before this step (6 x nIP)
    %   u: the displacement at all dofs before this step (2*ndof x 1)
    %   Ce: elastic tangent operator
    %   ALPHA: back stress for all IPs before this step (6 x nIP)
    %   R: radius for all IPs before this step (1 x nIP)
    %   mu: shear modulus
    %   H: hardening modulus
    %   beta: kinematic-isotropic interpolation factor (1=isotropic)
    %   f_ext: the vector of the applied load (2*ndof x 1)
    %   nodes: the nodal coordinates (nnode x 2)
    %   elements: the node numbers for each element (nel x 4)
    %   t: the plate thickness
    % Returns:
    %   SIG: the stress at all IPs after this step (6 x nIP)
    %   u: the displacement at all dofs after this step (2*ndof x 1)
    %   ALPHA: back stress for all IPs after this step (6 x nIP)
    %   R: radius for all IPs after this step (1 x nIP)

    u_diff=zeros(size(u,1),1);
    [f_int, A, SIG, ALPHA, R]=F_INT(SIG,u_diff,Ce,ALPHA,R,mu,H,beta, nodes, elements,t,nu);
    res=f_ext-f_int;

    % threshold for convergence
    tol=1e-9;

    SIG0=SIG;
    ALPHA0=ALPHA;
    R0=R;
    u_diff=zeros(size(u,1),1);

    ct=0;
    limit=200;
    while norm(res(free_dofs)) > tol && ct < limit
        ct = ct + 1;
        if (ct==limit)
            fprintf("Not converged: Newton iteration %d, residual = %.3e\n", ct, norm(res(free_dofs)));
        end
    
        % Update total increment
        u_diff(free_dofs) = u_diff(free_dofs) + A(free_dofs, free_dofs) \ res(free_dofs);
    
        % Evaluate new internal force and tangent
        [f_int, A, SIG, ALPHA, R] = F_INT(SIG0, u_diff, Ce, ALPHA0, R0, mu, H, beta, nodes, elements, t, nu);
        res = f_ext - f_int;
    end
    fprintf("Newton converged in %d iterations.\n", ct);
    u=u+u_diff;
end