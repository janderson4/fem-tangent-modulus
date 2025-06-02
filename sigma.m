function [SIG, ALPHA, R, C_alg] = sigma(SIG,Ce,EPSI,ALPHA,R,mu,H,beta)
    % Function to integrate stresses with isotropic and kinematic
    % hardening for a single point. 
    % Args:
    %   SIG: the stress at this point before this step (6 x nIP)
    %   Ce: elastic tangent operator
    %   EPSI: the change in strain at this step (6 x 1)
    %   ALPHA: back stress at this point before this step (6 x nIP)
    %   R: radius at this point before this step (1 x nIP)
    %   mu: shear modulus
    %   H: hardening modulus
    %   beta: kinematic-isotropic interpolation factor (1=isotropic)
    % Returns:
    %   SIG: the stress at this point after this step (6 x nIP)
    %   ALPHA: back stress at this point after this step (6 x nIP)
    %   R: radius at this point after this step (1 x nIP)

    % convert to tensorial shear strain
    EPSI(4:6)=EPSI(4:6)./2;

    % identity
    ID=[1;1;1;0;0;0];
    
    sig_tr=SIG+Ce*EPSI;
    
    p=1/3*sum(sig_tr(1:3));

    % CALCULATE TRIAL STRESSES

    s_tr=sig_tr-p*ID;
    XI=s_tr-ALPHA;
    
    % since the shear strain is engineering, it must be halved to get the
    % tensorial shear strain before squaring
    a=sqrt(sum(XI(1:3).^2)+2*sum((XI(4:6)).^2));
    
    % CHECK IF ELASTIC
    if (a <= R)
        SIG=sig_tr;
        C_alg=Ce;
        return
    end
    
    % PLASTIC PHASE: CALCULATE N (ALSO STORED IN "XI")
    nhat=XI/a;
    delta_lambda=(a-R)/(2*mu+H);
    
    SIG=sig_tr-2*mu*delta_lambda*nhat;
    
    R=R+beta*H*delta_lambda;
    
    ALPHA=ALPHA+(1-beta)*H*delta_lambda*nhat;
    
    Cep=Ce-2*mu*nhat*nhat'/(1+H/(2*mu));
    
    psi=2*mu*delta_lambda/a;
    
    Ia=eye(6);
    Ia(4:6,4:6)=Ia(4:6,4:6)-0.5*eye(3);

    C_alg=Cep-2*mu*psi*(Ia-1/3*ID*ID'-nhat*nhat');
end