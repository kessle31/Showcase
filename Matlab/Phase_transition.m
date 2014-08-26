%% Parameters


clear all
a       = 1;
gamma   = 0.1*a;
omega_0 = a/2;
Omega_0 = a/2;

%domega  = -0*omega_0; %delta omega
%omega   = omega_0 + domega; %hyperfine detuning

steps   = 50;
b       = [ 0 ; 0 ; -gamma/2];

syms X Y Z Sz Sx Sy bx by real;

%% Loop

for ii=1:steps
    
    Omega(ii)= 2 * Omega_0 * ii/steps + 0.0001;
    
    for jj=1:steps
    
        omega(jj) = 2 * omega_0 * jj/steps ;
        
    %% Semiclassical solution
    
    A=solve(-gamma*(Sz + 1/2) + 2*Omega(ii)*Sy,-gamma/2*Sx - omega(jj)*Y,...
        a*(Sy*Z-Sz*Y) - omega(jj)*Y,-gamma/2*Sy - 2*Omega(ii)*Sz + omega(jj)*X,...
        -a*(Sx*Z - Sz*X) + omega(jj)*X,Sy*X-Sx*Y,Z^2 + X^2 + Y^2 - 1);
    
    ASx=double(A.Sx);
    ASy=double(A.Sy);
    ASz=double(A.Sz);
    AX=double(A.X);
    AY=double(A.Y);
    AZ=double(A.Z);
    
%    Xi(ii,:)=AX;
%    Yi(ii,:)=AY;
%    Zi(ii,:)=AZ;
    
    
    if(size(ASx,1)~=2)
        warning('Number of solutions is not 2!')
    end
    
    %% Retrieving beta from the semicalssical solution
    
    B1=solve(bx*sqrt(2-bx^2-by^2) - real(AX(1)) , by*sqrt(2-bx^2-by^2) - (real(AY(1))-eps));
    B2=solve(bx*sqrt(2-bx^2-by^2) - real(AX(2)) , by*sqrt(2-bx^2-by^2) - (real(AY(2))-eps));

    
    
    %check which solution of beta gives the right value for Z
    if(double(abs(AZ(1)-(B1.bx(1)^2+B1.by(1)^2-1)))<10^(-10))
        
        B1bx=double(B1.bx(1));
        B1by=double(B1.by(1));
        
    elseif(double(abs(AZ(1)-(B1.bx(2)^2+B1.by(2)^2-1)))<10^(-10))
        
        B1bx=double(B1.bx(2));
        B1by=double(B1.by(2));
        
    else
        warning('no unique beta')
    end
    
    if(double(abs(AZ(2)-(B2.bx(1)^2+B2.by(1)^2-1)))<10^(-10))
        
        B2bx=double(B2.bx(1));
        B2by=double(B2.by(1));
        
    elseif(double(abs(AZ(2)-(B2.bx(2)^2+B2.by(2)^2-1)))<10^(-10))
        
        B2bx=double(B2.bx(2));
        B2by=double(B2.by(2));
        
    else
        warning('no unique beta')
    end
    
    %% examine the second order

    beta1(ii,jj)=B1bx - 1i *B1by;
    beta2(ii,jj)=B2bx - 1i *B2by;
    
    b1=beta1(ii,jj);
    b2=beta2(ii,jj);

    k1=2-abs(b1)^2;
    k2=2-abs(b2)^2;
    
    D_omega1=a*(abs(b1)^2-1);
    D_omega2=a*(abs(b2)^2-1);

    Omega_t1= Omega(ii) + a/2 *sqrt(k1) *b1;
    Omega_t2= Omega(ii) + a/2 *sqrt(k2) *b2;
    
    %optical Bloch matrices

    M1=[-gamma/2 + 1i*D_omega1              0                -2i*Omega_t1';...
                 0                -gamma/2 - 1i*D_omega1      2i*Omega_t1 ;...
          -1i*Omega_t1              1i*Omega_t1'                 -gamma];

    M2=[-gamma/2 + 1i*D_omega2              0                -2i*Omega_t2';...
                 0                -gamma/2 - 1i*D_omega2      2i*Omega_t2 ;...
            -1i*Omega_t2              1i*Omega_t2'                 -gamma];
        
    MI1=inv(M1);
    MI2=inv(M2);
    
    %steady state values of single electron operators Sp Sm Sz    
    vS1=-MI1*b;
    vS2=-MI2*b;
    
    %coefficients of operators A
        
    v1=[a/(4*sqrt(k1)) * (2*k1 - abs(b1)^2) ; -a/(4*sqrt(k1))*b1'^2 ; a*b1' ];
    vd1=[-a/(4*sqrt(k1))*b1^2 ; a/(4*sqrt(k1)) * (2*k1 - abs(b1)^2) ; a*b1 ];
   
    v2=[a/(4*sqrt(k2)) * (2*k2 - abs(b2)^2) ; -a/(4*sqrt(k2))*b2'^2 ; a*b2' ];
    vd2=[-a/(4*sqrt(k2))*b2^2 ; a/(4*sqrt(k2)) * (2*k2 - abs(b2)^2) ; a*b2 ];
    
    %steady state values of electron second moments
    
    
    mS1 = [ -vS1(1)^2                    vS1(3)+1/2-vS1(1)*vS1(2)      -vS1(1)*(1/2 + vS1(3));...
         1/2-vS1(3)-vS1(1)*vS1(2)            -vS1(2)^2                   vS1(2)*(1/2-vS1(3)) ;...
            vS1(1)*(1/2-vS1(3))             -vS1(2)*(1/2 + vS1(3))          1/4-vS1(3)^2   ];

    mS2 = [ -vS2(1)^2                    vS2(3)+1/2-vS2(1)*vS2(2)      -vS2(1)*(1/2 + vS2(3));...
         1/2-vS2(3)-vS2(1)*vS2(2)            -vS2(2)^2                   vS2(2)*(1/2-vS2(3)) ;...
            vS2(1)*(1/2-vS2(3))             -vS2(2)*(1/2 + vS2(3))          1/4-vS2(3)^2   ];  
        
    %coefficients
    
    F1a=-MI1*mS1;
    F1b=-mS1*MI1.';
    F2a=-MI2*mS2;
    F2b=-mS2*MI2.'; 
    
    AtA1   = v1.'*F1a*v1;
    AAt1   = v1.'*F1b*v1;
    AtdA1  = vd1.'*F1a*v1;
    AtAd1  = v1.'*F1a*vd1;
    
    AtA2   = v2.'*F2a*v2;
    AAt2   = v2.'*F2b*v2;
    AtdA2  = vd2.'*F2a*v2;
    AtAd2  = v2.'*F2a*vd2;
    
    Ra1=real(AtdA1);
    Rb1=real(AtAd1);
    Ia1=imag(AtdA1);
    Ib1=imag(AtAd1);
    c1=AtA1+AAt1;
    alpha1=(AtA1-AAt1)/2i;
    BB1=-(vS1(2)*(4*k1 + abs(b1)^2) + b1^2 * vS1(1)) * a*b1 /(16*sqrt(k1^3));
    F1 =-a*(4*k1 + abs(b1)^2)*(b1*vS1(1) + b1'*vS1(2))/(8*sqrt(k1^3)) + a*(vS1(3) + omega(jj)/a);
    
    Ra2=real(AtdA2);
    Rb2=real(AtAd2);
    Ia2=imag(AtdA2);
    Ib2=imag(AtAd2);
    c2=AtA2+AAt2;
    alpha2=(AtA2-AAt2)/2i;
    BB2=-(vS2(2)*(4*k2 + abs(b2)^2) + b2^2 * vS2(1)) * a*b2 /(16*sqrt(k2^3));
    F2 =-a*(4*k2 + abs(b2)^2)*(b2*vS2(1) + b2'*vS2(2))/(8*sqrt(k2^3)) + a*(vS2(3) + omega(jj)/a);
    
    G1=[ -((Ra1-Rb1) + 1i *(Ia1 +Ib1 + F1))         -2i*(alpha1' + BB1) ;...
               2i*(alpha1 + BB1')                   -((Ra1-Rb1) - 1i *(Ia1 +Ib1 + F1))  ];
 
    G2=[ -((Ra2-Rb2) + 1i *(Ia2 +Ib2 + F2))         -2i*(alpha2' + BB2) ;...
               2i*(alpha2 + BB2')                   -((Ra2-Rb2) - 1i *(Ia2 +Ib2 + F2))  ];      
           
    EV1(ii,jj,:)=eig(G1);
    EV2(ii,jj,:)=eig(G2);
    
    % Pic out the stable solution according to first moments
    if (EV1(ii,jj,1)<0 && EV1(ii,jj,2)<0)
        
        solswitch(ii,jj)=1;%solswitch pics out the stable semiclassical solution from A (cf. line 23)
        
        beta(ii,jj)=beta1(ii,jj);
        EV(ii,jj,:)=EV1(ii,jj,:);

        GG=[ -2*(Ra1-Rb1)+2i*(Ia1 +Ib1 + F1)                   0                 4i*(alpha1 + BB1') ;...
                     0                         -2*(Ra1-Rb1)-2i*(Ia1 +Ib1 + F1)    -4i*(alpha1' + BB1) ;...
                -2i*(alpha1' + BB1)                 2i*(alpha1 + BB1')                 -2*(Ra1-Rb1)];
    
        bb = [ -c1+2i*(alpha1+BB1') ; -c1'-2i*(alpha1'+BB1) ; 2*Rb1];

    elseif (EV2(ii,jj,1)<0 && EV2(ii,jj,2)<0)  
        
        solswitch(ii,jj)=2;
        
        beta(ii,jj)=beta2(ii,jj);       
        EV(ii,jj,:)=EV2(ii,jj,:);

      
        GG=[ -2*(Ra2-Rb2)+2i*(Ia2 +Ib2 + F2)                   0                 4i*(alpha2 + BB2') ;...
                        0                         -2*(Ra2-Rb2)-2i*(Ia2 +Ib2 + F2)    -4i*(alpha2' + BB2) ;...
                -2i*(alpha2' + BB2)                 2i*(alpha2 + BB2')                 -2*(Ra2-Rb2)]; 
      
        bb = [ -c2+2i*(alpha2+BB2') ; -c2'-2i*(alpha2'+BB2) ; 2*Rb2];  
    
        
    else     
        warning('no stable solution from first moments')
    
    end
    
    %steady state fluctuations
    fluct(ii,jj,:) = -GG\bb;
    
    covmat=[(fluct(ii,jj,1) + fluct(ii,jj,2) + 2*fluct(ii,jj,3) +1)            -1i*(fluct(ii,jj,2) - fluct(ii,jj,1)) ;...
                 -1i*(fluct(ii,jj,2) - fluct(ii,jj,1))              -(fluct(ii,jj,1) + fluct(ii,jj,2) - 2*fluct(ii,jj,3) - 1)];
    
    D(ii,jj)=sqrt(max(eig(covmat))*min(eig(covmat)));
    
    lammin(ii,jj)=real(min(eig(covmat)));
    
    %Semiclassical Observables 
    
    Xi(ii,jj)=AX(solswitch);
    Yi(ii,jj)=AY(solswitch);
    Zi(ii,jj)=AZ(solswitch);
    Xs(ii,jj)=ASX(solswitch);
    Ys(ii,jj)=ASY(solswitch);
    Zs(ii,jj)=ASZ(solswitch);    
    
    end
    
end