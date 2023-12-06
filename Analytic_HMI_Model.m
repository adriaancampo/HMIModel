% -------------------------------------------------------------------------
% Introduction: 
% -------------------------------------------------------------------------
% Disease processes often induce microstructural alterations in tissues,
% resulting in changes to the mechanical properties of the corresponding
% biological structures. Harmonic Motion Imaging (HMI) emerges as a 
% pivotal elasticity imaging technique adept at investigating the mechanical 
% parameters of tissues. By detecting tissue responses to an oscillatory 
% acoustic radiation force, HMI enables exploration of displacement, strain, 
% and Shear Wave Velocity (SWV). This code demonstrates an analytical HMI 
% model, facilitating swift and precise computations of these mechanical 
% parameters within a homogenous linear elastic material. The model offers 
% valuable insights into underlying biomechanical dynamics, with applications 
% extending to diverse scenarios, including tumor detection, characterization,
% and the monitoring of ablation procedures.
% -------------------------------------------------------------------------
% Usage:
% -------------------------------------------------------------------------
% - Modify parameters and settings as needed.
% - Run the script to execute the simulations.
% -------------------------------------------------------------------------
% Output:
% -------------------------------------------------------------------------
% - The code generates displacement and strain field components as matrices,
%   representing the simulated spatial distribution.
% -------------------------------------------------------------------------
% Authors: 
% -------------------------------------------------------------------------
% Dr. Adriaan Campo, Dr. Matt McGarry
% Date: 06/12/2023
% -------------------------------------------------------------------------
% Cite as: 
% -------------------------------------------------------------------------
% Matthew D J McGarry et al 2021 Phys. Med. Biol. 66 075017 
% 
% An analytical model of full-field displacement and strain induced by 
% amplitude-modulated focused ultrasound in harmonic motion imaging
% Matthew D J McGarry, Adriaan Campo
% Published 6 April 2021
% Physics in Medicine & Biology, Volume 66, Number 7 
% DOI 10.1088/1361-6560/abddd1  
% -------------------------------------------------------------------------

%% Prepare workspace
clc
clear
close all

%% Define parameters and settings
res_force_volume = 101; % resolution force volume
res_total_volume = 61; % resolution total volume
FOV = 0.003; % radius total volume (m)
freq = 50; % driving frequency (1/s)
rho = 1000; % density matrix (kg/m^3)
v=0.49; % poisson ratio matrix 
E = 10000; % young's modulus matrix (N/m^2)
mu=E./(2.*(1+v)); % mu value matrix (N/m^2)
M = E.*(1-v)./((1+v).*(1-2.*v)); % M value matrix (N/m^2)
omg=2.*pi.*freq; % radial frequency (rad/s)
kc=sqrt(rho.*(2.*pi.*freq).^2./M); % kc value matrix (kg/(N*m^2*s^2))^(1/2)
ks=sqrt(rho.*(2.*pi.*freq).^2./mu); % ks value matrix (kg/(N*m^2*s^2))^(1/2)
BF = -10000; % force amplitude 

%% Setup main grid for calculations
xp_force_volume = linspace(-FOV,FOV,res_force_volume);
yp_force_volume = linspace(-FOV,FOV,res_force_volume);
zp_force_volume = linspace(-FOV,FOV,res_force_volume);

xp_total_volume = linspace(-FOV,FOV,res_total_volume);
yp_total_volume = linspace(-FOV,FOV,res_total_volume);
zp_total_volume = linspace(-FOV,FOV,res_total_volume);

S1 = find(xp_total_volume==0);
S2 = find(yp_total_volume==0);

X_force_volume = xp_force_volume;
Y_force_volume = zp_force_volume;
Z_force_volume = yp_force_volume;
[XX_force_volume,YY_force_volume,ZZ_force_volume] = ndgrid(X_force_volume,Y_force_volume,Z_force_volume);

X_total_volume = xp_total_volume;
Y_total_volume = zp_total_volume;
Z_total_volume = yp_total_volume;
[XX_total_volume,YY_total_volume,ZZ_total_volume] = ndgrid(X_total_volume,Y_total_volume,Z_total_volume);

%% Setup subgrid for Singularity Avoidance
dL = (xp_force_volume(2)-xp_force_volume(1));
dF = BF.*dL.*dL.*dL;
subgrid_z = [dL./4, dL./4,dL./4, dL./4,-dL./4,-dL./4,-dL./4,-dL./4];
subgrid_x = [dL./4,-dL./4,dL./4,-dL./4, dL./4,-dL./4, dL./4,-dL./4];
subgrid_y = [dL./4,-dL./4,-dL./4,dL./4, dL./4,-dL./4,-dL./4, dL./4];
% % above subgrid definition is better
% subgrid_z = [dL./4,-dL./4];
% subgrid_x = [0,0];
% subgrid_y = [0,0];

%% Define force volume cylinder Coordinates         
Zco = find(ZZ_force_volume<=0.0015&ZZ_force_volume>=-0.0015);
Rgrid = sqrt(XX_force_volume.^2 + YY_force_volume.^2);
[ForceVolumeCoos, ~]=find(Rgrid(:)<0.00075);
[Xco, Yco]=ind2sub(size(Rgrid),ForceVolumeCoos);
ForceVolumeCoos = intersect(Zco,ForceVolumeCoos);

%% Simulation Results Initialization
S = round(size(XX_total_volume,2)/2);
X_FOV = XX_total_volume(S2:end,S1,:);
Y_FOV = YY_total_volume(S2:end,S1,:);
Z_FOV = ZZ_total_volume(S2:end,S1,:);

dFsub = dF./length(subgrid_z);
cte1 = (dFsub).*1./(4.*pi.*rho.*omg.^2);

resultUx = 0.*X_FOV;
resultUy = 0.*X_FOV;
resultUz = 0.*X_FOV;

resultUz_near_S = 0.*X_FOV;
resultUz_intermediate_S = 0.*X_FOV;
resultUz_far_S = 0.*X_FOV;

resultUz_near_P = 0.*X_FOV;
resultUz_intermediate_P = 0.*X_FOV;
resultUz_far_P = 0.*X_FOV;

resultEx = 0.*X_FOV;
resultEy = 0.*X_FOV;
resultEz = 0.*X_FOV;
resultExz = 0.*X_FOV;
resultEyz = 0.*X_FOV;

%% Main Simulation Loop
for idx1 = 1:length(ForceVolumeCoos)
    
    %% Loop through force volume coordinates
    Xcenter = XX_force_volume(ForceVolumeCoos(idx1));
    Ycenter = YY_force_volume(ForceVolumeCoos(idx1));
    Zcenter = ZZ_force_volume(ForceVolumeCoos(idx1));
    
    for idx2 = 1:length(subgrid_z)

        %% Loop through subgrid coordinates
        Xsubcenter = Xcenter + subgrid_x(idx2);
        Ysubcenter = Ycenter + subgrid_y(idx2);
        Zsubcenter = Zcenter + subgrid_z(idx2);

        xk=(X_FOV-Xsubcenter);
        yk=(Y_FOV-Ysubcenter);
        zk=(Z_FOV-Zsubcenter);

        %% Ux and Uy field
        rk = sqrt(xk.^2 + yk.^2 + zk.^2);
        zr=zk.^2./rk.^2;

        C1 = ( kc.^2.*zk./(rk.^3) + kc.*(3.*1i.*zk)./(rk.^4) - (3.*zk)./(rk.^5) ) .* exp(1i.*kc.*rk);
        C2 = (-ks.^2.*zk./(rk.^3) - ks.*(3.*1i.*zk)./(rk.^4) + (3.*zk)./(rk.^5) ) .* exp(1i.*ks.*rk);

        uxt = cte1 .* xk .* C1 ...
            + cte1 .* xk .* C2 ;

        uyt = cte1 .* yk .* C1 ...
            + cte1 .* yk .* C2 ;

        %% Uz S-field
        uzt_near_s         = cte1 .* ( (3.*zr-1)./(rk.^3)          ) .* exp(1i.*ks.*rk);                                     
        uzt_intermediate_s = cte1 .* ( ks.*(1i-3.*1i.*zr)./(rk.^2) ) .* exp(1i.*ks.*rk);
        uzt_far_s          = cte1 .* ( ks.^2.*(1-zr)./rk           ) .* exp(1i.*ks.*rk);

        %% Uz P-field
        uzt_near_p         = cte1 .* ( (1-3.*zr)./(rk.^3)          ) .* exp(1i.*kc.*rk);
        uzt_intermediate_p = cte1 .* ( kc.*(3.*1i.*zr-1i)./(rk.^2) ) .* exp(1i.*kc.*rk);
        uzt_far_p          = cte1 .* ( kc.^2.*zr./rk               ) .* exp(1i.*kc.*rk);

        %% Uz total field
        uzt = uzt_far_p ...
            + uzt_intermediate_p ...
            + uzt_near_p ...
            ...
            + uzt_far_s ...
            + uzt_intermediate_s ...
            + uzt_near_s;                      

        %% Sum fields for total solution
        C3 = (- (cte1.*zk.*exp(kc.*rk.*1i).*(3.*kc.^2.*rk.^2 + kc.*rk.*12i - 15))./(rk.^7) + (cte1.*zk.*exp(ks.*rk.*1i).*(3.*ks.^2.*rk.^2 + ks.*rk.*12i - 15))./(rk.^7) + (cte1.*kc.*zk.*exp(kc.*rk.*1i).*(kc.^2.*rk.^2 + kc.*rk.*3i - 3).*1i)./(rk.^6) - (cte1.*ks.*zk.*exp(ks.*rk.*1i).*(ks.^2.*rk.^2 + ks.*rk.*3i - 3).*1i)./(rk.^6));
        C4 = (cte1.*zk.*exp(kc.*rk.*1i).*(kc.^2.*rk.^4 + kc.*rk.^3.*3i - 3.*rk.^2))./(rk.^7) - (cte1.*zk.*exp(ks.*rk.*1i).*(ks.^2.*rk.^4 + ks.*rk.^3.*3i - 3.*rk.^2))./(rk.^7);
        C5 = (2.*cte1.*(kc.*exp(kc.*rk.*1i).*3i - ks.*exp(ks.*rk.*1i).*3i + kc.^3.*zk.^2.*exp(kc.*rk.*1i).*1i - ks.^3.*zk.^2.*exp(ks.*rk.*1i).*1i))./rk.^4 + (cte1.*(2.*kc.^2.*exp(kc.*rk.*1i) - 3.*ks.^2.*exp(ks.*rk.*1i)))./rk.^3 - (6.*cte1.*(exp(kc.*rk.*1i) - exp(ks.*rk.*1i) + 2.*kc.^2.*zk.^2.*exp(kc.*rk.*1i) - 2.*ks.^2.*zk.^2.*exp(ks.*rk.*1i)))./rk.^5 + (cte1.*ks.^3.*exp(ks.*rk.*1i).*1i)./rk.^2 + (30.*cte1.*zk.^2.*(exp(kc.*rk.*1i) - exp(ks.*rk.*1i)))./rk.^7 - (30.*cte1.*zk.^2.*(kc.*exp(kc.*rk.*1i).*1i - ks.*exp(ks.*rk.*1i).*1i))./rk.^6;

        ext = C3.*xk.^2 + C4;
        eyt = C3.*yk.^2 + C4;
        ezt = (cte1.*zk.*exp(kc.*rk.*1i).*(2.*kc.^2.*rk.^4 - 3.*kc.^2.*rk.^2.*zk.^2 + kc.*rk.^3.*8i - kc.*rk.*zk.^2.*12i - 9.*rk.^2 + 15.*zk.^2))./(rk.^7) ...
            - (cte1.*zk.*exp(ks.*rk.*1i).*(3.*ks.^2.*rk.^4 - 3.*ks.^2.*rk.^2.*zk.^2 + ks.*rk.^3.*8i - ks.*rk.*zk.^2.*12i - 9.*rk.^2 + 15.*zk.^2))./(rk.^7) ...
            + (cte1.*kc.*zk.*exp(kc.*rk.*1i).*(kc.^2.*rk.^2.*zk.^2 - kc.*rk.^3.*1i + kc.*rk.*zk.^2.*3i + rk.^2 - 3.*zk.^2).*1i)./(rk.^6) ...
            + (cte1.*ks.*zk.*exp(ks.*rk.*1i).*(ks.^2.*rk.^4 - ks.^2.*rk.^2.*zk.^2 + ks.*rk.^3.*1i - ks.*rk.*zk.^2.*3i - rk.^2 + 3.*zk.^2).*1i)./(rk.^6);

        exzt = C5.*xk;
        eyzt = C5.*yk;

        resultUx = resultUx + uxt;
        resultUy = resultUy + uyt;
        resultUz = resultUz + uzt;
        
        resultUz_near_S = resultUz_near_S + uzt_near_s;
        resultUz_intermediate_S = resultUz_intermediate_S + uzt_intermediate_s;
        resultUz_far_S = resultUz_far_S + uzt_far_s;
        
        resultUz_near_P = resultUz_near_P + uzt_near_p;
        resultUz_intermediate_P = resultUz_intermediate_P + uzt_intermediate_p;
        resultUz_far_P = resultUz_far_P + uzt_far_p;

        resultEx = resultEx + ext;
        resultEy = resultEy + eyt;
        resultEz = resultEz + ezt;
        resultExz = resultExz + exzt;
        resultEyz = resultEyz + eyzt;

    end 

end

%% Plot figures
%% Plot simulation grid
figure
hold on
plot3(XX_force_volume(1:37:end),YY_force_volume(1:37:end),ZZ_force_volume(1:37:end),'.','Linewidth',0.1)
plot3(XX_force_volume(ForceVolumeCoos),YY_force_volume(ForceVolumeCoos),ZZ_force_volume(ForceVolumeCoos),'ro','Linewidth',2)
title('simulation grid')
xlabel('x position (m)')
ylabel('y position (m)')
zlabel('z position (m)')
legend({'overall grid','force volume grid'})
view([45 45])
axis square

%% Plot displacement fields
figure
subplot(3,1,1)
imagesc(squeeze(real(resultUx))')
axis square
title('U_x displacement')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,1,2)
imagesc(squeeze(real(resultUy))')
axis square
title('U_y displacement')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,1,3)
imagesc(squeeze(real(resultUz))')
axis square
title('U_z displacement')
xlabel('x position (m)')
ylabel('z position (m)')

%% Plot strain fields
figure
subplot(3,2,1)
imagesc(squeeze(real(resultEx))')
axis square
title('E_x strain')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,3)
imagesc(squeeze(real(resultEy))')
axis square
title('E_y strain')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,5)
imagesc(squeeze(real(resultEz))')
axis square
title('E_z strain')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,2)
imagesc(squeeze(real(resultExz))')
axis square
title('E_x_z strain')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,4)
imagesc(squeeze(real(resultEyz))')
axis square
title('E_y_z strain')
xlabel('x position (m)')
ylabel('z position (m)')

%% Plot S and P displacement fields, for near, intermediate and far fields
figure
subplot(3,2,1)
imagesc(squeeze(real(resultUz_near_S))')
axis square
title('U_z near S')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,3)
imagesc(squeeze(real(resultUz_intermediate_S))')
axis square
title('U_z intermediate S')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,5)
imagesc(squeeze(real(resultUz_far_S))')
axis square
title('U_z far S')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,2)
imagesc(squeeze(real(resultUz_near_P))')
axis square
title('U_z near P')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,4)
imagesc(squeeze(real(resultUz_intermediate_P))')
axis square
title('U_z intermediate P')
% xlabel('x position (m)')
ylabel('z position (m)')
subplot(3,2,6)
imagesc(squeeze(real(resultUz_far_P))')
axis square
title('U_z far P')
xlabel('x position (m)')
ylabel('z position (m)')
