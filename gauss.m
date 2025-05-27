function [xi, eta] = gauss(index)
    % Gauss points and weights for 2x2 quadrature
    % Args:
    %   index: the number of the Gauss point (1-4)
    % Returns:
    %   xi: the xi coordinate
    %   eta: the eta coordinate

    if (mod(index,2)==1)
        xi=-1/sqrt(3);
    else
        xi=1/sqrt(3);
    end
    
    if (index<3)
        eta=-1/sqrt(3);
    else
        eta=1/sqrt(3);
    end
end