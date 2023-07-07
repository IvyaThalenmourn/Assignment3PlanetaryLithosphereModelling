%% Assignment 3
% Class: AE4893 - Physics of Planetary Interiors
% Author: Alyona Glazyrina
% Date: 3 July 2023
% Requires: Mapping Toolbox, GSH code

%axesm('MapProjection','mollweid', 'Frame', 'on', 'Grid', 'on');
%geoshow(flip(t_crust)./1e3, [1 90 180], 'DisplayType', 'texturemap');

addpath('Tools/')
% Read topography data
resolution = 4;
f = fopen('megt90n000cb.img','r','ieee-be');
topo = fread(f,[360*resolution Inf],'int16')';
topo = imresize(topo,[180,360]);
fclose(f); 

% Read gravity data
load MarsGravity.txt
V = [0 0 1 0; MarsGravity(:,1:4)];
V = sortrows(V,2); V(1,3)=0; V(3,3)=0;
Model = struct();
Model.Re = 3396000;
Model.Re_analyse = 3396000;
Model.GM = 42828.3748574*10^9;
grav_data = model_SH_synthesis([0.5 359.5 1],[-89.5 89.5 1],0,[0 120],V,Model);

% Miscellaneous
Model.nmax = 120;
Model.number_of_layers = 2;
Model.geoid = 'none';
Model.rho_c = 2700; % [kg/m^3] - source: 
Model.rho_m = 3550; % [kg/m^3] - source: 
D = 50000; % [m] - reference thickness - source: 

clear MarsGravity f ans resolution

%% Question 2
% Topography Visualization
lon = grav_data.grd.lon(1,:);
lat = grav_data.grd.lat(:,1);

figure(1)
subplot(2,2,1)
imagesc(lon,flip(lat),topo./1e3);c=colorbar; title("Topography");
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Topography (km)'); set(gca,'YDir','normal'); colormap(jet);
yticks([-60 -30 0 30 60]); xticks([-90 0 90])

% Gravity Visualization
grav_obs = gradient(flip(grav_data.pot))/(Model.GM/Model.Re^2);

subplot(2,2,2)
imagesc(lon,lat,flipud(grav_obs)); c=colorbar; 
xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)]);
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
yticks([-60 -30 0 30 60]); xticks([-90 0 90]);
title('Potential gravity field')
ylabel(c,'mGal'); set(gca,'YDir','normal');

% Bouguer Anomaly
subplot(2,2,3)
dataG = flip(grav_data.vec.R);
imagesc(lon,lat,flip(dataG)*1e5); c=colorbar; 
xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)])
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)')
yticks([-60 -30 0 30 60]); xticks([-90 0 90])
title('R-vector w/out Bouguer Correction')
ylabel(c,'mGal'); set(gca,'YDir','normal');

% Bouguer Correction
subplot(2,2,4)
G = 6.67430*10^-11; % [N*m^2/kg^2]
grav_bouguer = flip(grav_data.vec.R)-2*pi*Model.rho_c*G*topo;
imagesc(lon,lat,flip(grav_bouguer*1e5));c=colorbar; 
xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)])
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)')
yticks([-60 -30 0 30 60]); xticks([-90 0 90])
title('R-vector w/ Bouguer Correction')
ylabel(c,'mGal'); set(gca,'YDir','normal');

%% Question 3 - Bouguer Inversion
dataG = imresize(dataG,[45,90]);
dataG = imresize(dataG,[180,360]);

Model.moho = -D*ones(size(topo));

k = 15;

err = 1;
iter = 1;
threshold = 0.1;
residual = 0;
while err>threshold && iter<11

    Model.moho = Model.moho+k*residual*10^5;
    for i = 1:size(Model.moho,1)
        for j = 1:size(Model.moho,2)
            if Model.moho(i,j) > topo(i,j)
                Model.moho(i,j) = topo(i,j);
            end
        end
    end

    V_boug = segment_2layer_model(topo,Model.moho,-500E3,Model.rho_c,Model.rho_m,25E3,Model);
    V_boug(1,3) = 0; V_boug(2,2) = 0; V_boug(3,3) = 0;

    [model_data] = model_SH_synthesis([0.5 359.5 1],[-89.5 89.5 1],0,[0 50],V_boug,Model);

    residual = (dataG - flip(model_data.vec.R));
    if iter == 1
        residual0 = residual;
    end

    %err=abs(mean(residual,"all"))*10^5;
    err = max(max(abs(residual)))*10^5;
    fprintf("Iteration %.0f error: %.2f.\n",iter,err)

    iter = iter+1;
end

%%
figure(3)
subplot(3,2,1)
imagesc(lon,flip(lat),dataG*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('DataG'); set(gca,'YDir','normal')
subplot(3,2,2)
imagesc(lon,flip(lat),flip(model_data.vec.R)*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('ModelG'); set(gca,'YDir','normal')
subplot(3,2,3)
imagesc(lon,flip(lat),topo/1e3); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'km'); title('Topography'); set(gca,'YDir','normal')
subplot(3,2,4)
imagesc(lon,flip(lat),Model.moho/1e3); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'km'); title('Crust-mantle boundary'); set(gca,'YDir','normal')
subplot(3,2,5)
imagesc(lon,flip(lat),residual0*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('Initial residual');
set(gca,'YDir','normal')
subplot(3,2,6)
imagesc(lon,flip(lat),residual*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('Residual after x iterations')
set(gca,'YDir','normal')

%% Question 4 - Airy Model
r_airy = Model.rho_c*topo/(Model.rho_m-Model.rho_c); % Depth update
t_crust = topo+D*ones(size(topo))+r_airy;
airy_mb = -D*ones(size(topo))-r_airy;

airy_sc = cs2sc(GSHA(t_crust,120));

figure(3)
imagesc(lon,flip(lat),-t_crust./1e3); c=colorbar; colormap(jet)
xlabel('Longitude (\circ)');ylabel('Latitude (\circ)');title('Airy model');
ylabel(c,'Crustal depth (km)'); set(gca,'YDir','normal');

V_airy = segment_2layer_model(topo,airy_mb,-300E3,Model.rho_c,Model.rho_m,25E3,Model);

clear r_airy t_crust

%% Question 5 - Flexural Model
Te = 95*10^3; % [m] - source: 
[n,phi,phi2,flex_mb2] = flexural_model(Te,Model,airy_sc,lon,lat);

figure(4)
semilogx(n,phi); hold on; semilogx(n,phi2)
xlim([2 180]); xticks([2 5 10 25 50 100])
xlabel('Spherical Harmonic Degree (n)'); ylabel('Flexural Response')
title('Flexural response function'); legend("Flat plate","Shell")

figure(5)
imagesc(lon,flip(lat),-flex_mb2/10^3); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)')
ylabel(c,'Crustal depth (km)'); title('Flexural model (shell)')
set(gca,'YDir','normal')
colormap(jet);

V_flex = segment_2layer_model(topo,flex_mb2,-300E3,Model.rho_c,Model.rho_m,25E3,Model);

clear phi phi2 G Te

%% Question 6 - Gravity Data Insertion
[~,DV_0] = degreeVariance(V);
[~,DV_2] = degreeVariance(V_airy);
[~,DV_4] = degreeVariance(V_flex);

figure(6)
loglog(n,DV_0*10^10,'*','MarkerSize',4)
hold on
loglog(n,DV_2*10^10,'*','MarkerSize',4)
loglog(n,DV_4*10^10,'*','MarkerSize',4)
xlabel('Spherical Harmonic Degree (n)'); ylabel('Power spectrum');
xlim([2 120]); xticks([2 5 10 25 50 100]); grid on
legend('Measured','Airy','Flexural [shell]');

%% Question 7 - Optimal T_e
Te = linspace(40e3,150e3,12); % reference: 
rmse_array = zeros(1,length(Te));
DV_01 = DV_0(2:50);
for i = 1:length(Te)
    fprintf("Te = %.0f km.\n",Te(i)/1e3)
    [~,~,~,flex_mb] = flexural_model(Te(i),Model,airy_sc,lon,lat);
    V_flex = segment_2layer_model(topo,flex_mb,-300E3,Model.rho_c,Model.rho_m,25E3,Model);
    V_flex(1,3) = 0; V_flex(2,2) = 0; V_flex(3,3) = 0;
    [~,DV_5] = degreeVariance(V_flex);
    rmse_array(i) = rmse(DV_0(2:90),DV_5(2:90));
end

figure(7)
plot(Te/1e3,rmse_array*10^10);
xlabel('Te (km)'); ylabel('RMS error (mGal)');
title("Root mean square error for various Te values")
%%
[n,phi,phi2,flex_mb2] = flexural_model(120e3,Model,airy_sc,lon,lat);

figure(8)
semilogx(n,phi)
hold on
semilogx(n,phi2)
xlim([2 180]); xticks([2 5 10 25 50 100])
xlabel('Spherical Harmonic Degree (n)'); ylabel('Flexural Response Function')
title('Flexural response function'); legend("Flat plate","Shell")
colormap(jet);
figure(9)
imagesc(lon,flip(lat),-flex_mb2/10^3);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)')
ylabel(c,'Crustal depth (km)'); title('Flexural model (shell)')
set(gca,'YDir','normal')
colormap(jet);

%% Question 8 - Lateral Density Variations
Model.moho = -flex_mb2;

k = 15;

err = 1;
iter = 1;
threshold = 0.1;
residual = 0;
while err>threshold && iter<3

    Model.rho_c = Model.rho_c+k*residual*10^5;
    
    V = segment_2layer_model(topo,Model.moho,-500E3,Model.rho_c,Model.rho_m,25E3,Model);
    V(1,3) = 0; V(2,2) = 0; V(3,3) = 0;

    [model_data] = model_SH_synthesis([0.5 359.5 1],[-89.5 89.5 1],0,[0 50],V,Model);

    residual = (dataG - flip(model_data.vec.R));
    if iter == 1
        residual0 = residual;
    end

    err=abs(mean(residual,"all"))*10^5;
    fprintf("Iteration %.0f error: %.2f.\n",iter,err)

    iter = iter+1;
end

figure(10)
figure(3)
subplot(3,2,1)
imagesc(lon,flip(lat),dataG*10^5);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal');
title('DataG')
set(gca,'YDir','normal')
%colormap(jet)
subplot(3,2,2)
imagesc(lon,flip(lat),flip(model_data.vec.R)*10^5);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal');
title('ModelG')
set(gca,'YDir','normal')
%colormap(jet)
subplot(3,2,3)
imagesc(lon,flip(lat),topo/1e3);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)'); ylabel(c,'km');
title('Topo')
set(gca,'YDir','normal')
%colormap(jet)
subplot(3,2,4)
imagesc(lon,flip(lat),Model.rho_c/1e3);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)'); ylabel(c,'km');
title('Density')
set(gca,'YDir','normal')
%colormap(jet)
subplot(3,2,5)
imagesc(lon,flip(lat),residual0*10^5);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal');
title('Initial residual')
set(gca,'YDir','normal')
%colormap(jet)
subplot(3,2,6)
imagesc(lon,flip(lat),residual*10^5);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal');
title('Residual after x iterations')
set(gca,'YDir','normal')
%colormap(jet)

%% Flexural Model Calculation
function[n,phi,phi2,flex_mb2] = flexural_model(Te,Model,airy_sc,lon,lat)
    E = 65*10^9; % [Pa] - source: Crane & Rich, 20
    v = 0.25; % [poisson's ratio] - source: Crane & Rich, 20
    R = Model.Re;
    n = 1:size(airy_sc,1);
    D = E*Te^3/(12*(1-v^2));
    phi = (1 + D./((Model.rho_m-Model.rho_c)*3.721).*(2.*(n+1)./(2*R)).^4).^(-1);
    term2 = (12*(1-v^2)/Te^2/R^2).*(1-2./(n.*(n+1)))./(1-1./(n.*(n+1)));
    term1 = (n.*(n+1)-2).^2./(R^4.*(1-(1-v)./(n.*(n+1))));
    phi2 = (1+D./((Model.rho_m-Model.rho_c)*3.721).*(term1+term2)).^-1;
    
    sc_flex = zeros(size(airy_sc));
    sc_flex2 = zeros(size(airy_sc));
    
    for m = 1:size(airy_sc,2)
        sc_flex(:,m) = airy_sc(:,m).*phi';
        sc_flex2(:,m) = airy_sc(:,m).*phi2';
    end
    
    flex_mb2 = GSHS(sc_flex2,lon,90-flip(lat),120);
end

