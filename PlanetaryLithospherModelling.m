%% Assignment 3
% Class: AE4893 - Physics of Planetary Interiors
% Author: Alyona Glazyrina
% Date: 3 July 2023
% Requires: Mapping Toolbox, GSH code

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
Model.rho_c = 2700; % [kg/m^3] - source: Wieczorek & Zuber, 2004
Model.rho_m = 3550; % [kg/m^3] - source: Wieczorek & Zuber, 2004
D = 50000; % [m] - reference thickness - source: Wieczorek & Zuber, 2004

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

% Bouguer Anomaly Visualization
subplot(2,2,3)
imagesc(lon,lat,grav_data.vec.R*1e5); c=colorbar; 
xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)])
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)')
yticks([-60 -30 0 30 60]); xticks([-90 0 90])
title('R-vector w/out Bouguer Correction')
ylabel(c,'mGal'); set(gca,'YDir','normal');

% Bouguer Correction Visualization
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
% scaling down of gravity data
dataG = flip(grav_data.vec.R);
dataG = imresize(dataG,[45,90]);
dataG = imresize(dataG,[180,360]);

Model.rho_c = 2700; % [kg/m^3]
Model.moho = -D*ones(size(topo)); % crust-mantle boundary

k = 15; % scaling factor

err = 1000;
iter = 1;
threshold = 250; % [mGal]
residual = 0;
while err>threshold && iter<31

    Model.moho = Model.moho+k*residual*10^5; % moho adjustment
    % moho cannot exceed topo - unphysical
    for i = 1:size(Model.moho,1)
        for j = 1:size(Model.moho,2)
            if Model.moho(i,j) > topo(i,j)
                Model.moho(i,j) = topo(i,j);
            end
        end
    end

    % V of model
    V_boug = segment_2layer_model(topo,Model.moho,-500E3,Model.rho_c,Model.rho_m,25E3,Model);
    V_boug(1,3) = 0; V_boug(2,2) = 0; V_boug(3,3) = 0;

    % gravity data of model
    [model_data] = model_SH_synthesis([0.5 359.5 1],[-89.5 89.5 1],0,[0 50],V_boug,Model);

    % residual and error calcuation
    residual = (dataG - flip(model_data.vec.R));
    if iter == 1
        residual0 = residual;
    end
    err = max(max(abs(residual)))*10^5;
    fprintf("Iteration %.0f error: %.2f.\n",iter,err)

    iter = iter+1;
end

% Visualization
figure(3)
subplot(3,2,1) % observed gravity
imagesc(lon,flip(lat),dataG*10^5); c=colorbar; colormap(jet)
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('DataG'); set(gca,'YDir','normal')
subplot(3,2,2) % modeled gravity
imagesc(lon,flip(lat),flip(model_data.vec.R)*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('ModelG'); set(gca,'YDir','normal')
subplot(3,2,3) % topography
imagesc(lon,flip(lat),topo/1e3); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Topography (km)'); title('Topography'); set(gca,'YDir','normal')
subplot(3,2,4) % crustal depth
imagesc(lon,flip(lat),Model.moho/1e3); c=colorbar; set(gca,'YDir','normal')
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Crustal Depth (km)'); title('Crustal depth');
c.Ticks = [-140 -130 -120 -110 -100 -90 -80 -70 -60 -50 -40 -30 -20 -10];
subplot(3,2,5)
imagesc(lon,flip(lat),residual0*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('Initial residual');
set(gca,'YDir','normal')
subplot(3,2,6)
imagesc(lon,flip(lat),residual*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('Final residual')
set(gca,'YDir','normal');

clear G i j iter k threshold err residual residual0 model_data

%% Question 4 - Airy Model
r_airy = Model.rho_c*topo/(Model.rho_m-Model.rho_c); % Depth update
t_crust = topo+D*ones(size(topo))+r_airy;
airy_mb = -D*ones(size(topo))-r_airy;

airy_sc = cs2sc(GSHA(t_crust,120));

figure(3)
imagesc(lon,flip(lat),-t_crust./1e3); c=colorbar; colormap(jet)
xlabel('Longitude (\circ)');ylabel('Latitude (\circ)');title('Airy model');
ylabel(c,'Crustal depth (km)'); set(gca,'YDir','normal');
c.Ticks = [-140 -130 -120 -110 -100 -90 -80 -70 -60 -50 -40 -30 -20 -10];

V_airy = segment_2layer_model(topo,airy_mb,-500E3,Model.rho_c,Model.rho_m,25E3,Model);

clear r_airy

%% Question 5 - Flexural Model
Te = 150*10^3; % [m] - source: Thiriet, Michaut, Breuer, & Plesa, 2018
[n,phi,phi2,flex_mb2] = flexural_model(Te,Model,airy_sc,lon,lat);

figure(4) % spheric harmonic response visualization
semilogx(n,phi); hold on; semilogx(n,phi2)
xlim([2 180]); xticks([2 5 10 25 50 100])
xlabel('Spherical Harmonic Degree (n)'); ylabel('Flexural Response')
title('Flexural response function'); legend("Flat plate","Shell")

figure(5) % crustal depth visualization
imagesc(lon,flip(lat),-flex_mb2/10^3); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)')
ylabel(c,'Crustal depth (km)'); title('Flexural model (shell)')
set(gca,'YDir','normal')
colormap(jet);

V_flex = segment_2layer_model(topo,flex_mb2,-500E3,Model.rho_c,Model.rho_m,25E3,Model);

clear phi phi2 Te

%% Comparison with Literature
% Bouguer inversion
figure(6)
axesm('MapProjection','robinson', 'Frame', 'on', 'Grid', 'on');
temp_moho = -flipud(Model.moho); moho = zeros(size(Model.moho));
moho(:,1:180) = temp_moho(:,181:360); moho(:,181:360) = temp_moho(:,1:180);
geoshow(moho./1e3, [1 90 180], 'DisplayType', 'surface'); c=colorbar;
c.Ticks = [10 20 30 40 50 60 70 80]; ylabel(c,'Crustal depth (km)');
colormap(jet);
% Airy
figure(7)
axesm('MapProjection','robinson', 'Frame', 'on', 'Grid', 'on');
temp_moho = flipud(t_crust); moho = zeros(size(Model.moho));
moho(:,1:180) = temp_moho(:,181:360); moho(:,181:360) = temp_moho(:,1:180);
geoshow(moho./1e3, [1 90 180], 'DisplayType', 'surface'); c=colorbar;
c.Ticks = [10 20 30 40 50 60 70 80 90 100 110 120 130 140];
ylabel(c,'Crustal depth (km)');
colormap(jet);
% Flexural shell
figure(8)
axesm('MapProjection','robinson', 'Frame', 'on', 'Grid', 'on');
temp_moho = flipud(flex_mb2); moho = zeros(size(Model.moho));
moho(:,1:180) = temp_moho(:,181:360); moho(:,181:360) = temp_moho(:,1:180);
geoshow(moho./1e3, [1 90 180], 'DisplayType', 'surface'); c=colorbar;
c.Ticks = [10 20 30 40 50 60 70 80]; ylabel(c,'Crustal depth (km)');
colormap(jet);

%% Question 6 - Gravity Data Insertion
[~,DV_0] = degreeVariance(V);
[~,DV_1] = degreeVariance(V_boug);
[~,DV_2] = degreeVariance(V_airy);
[~,DV_4] = degreeVariance(V_flex);

figure(9)
loglog(n,DV_0*10^10,'*','MarkerSize',4)
hold on
loglog(n,DV_1*10^10,'*','MarkerSize',4)
loglog(n,DV_2*10^10,'*','MarkerSize',4)
loglog(n,DV_4*10^10,'*','MarkerSize',4)
xlabel('Spherical Harmonic Degree (n)'); ylabel('Power spectrum');
xlim([2 120]); xticks([2 5 10 25 50 100]); grid on
legend('Measured','Bouguer Inversion','Airy','Flexural [shell]');

%% Question 7 - Optimal T_e
Te = linspace(40e3,150e3,12); % reference: Thiriet, Michaut, Breuer, & Plesa, 2012

rmse_array = zeros(1,length(Te));
for i = 1:length(Te)
    fprintf("Te = %.0f km.\n",Te(i)/1e3)
    [~,~,~,flex_mb2] = flexural_model(Te(i),Model,airy_sc,lon,lat);
    V_flex = segment_2layer_model(topo,flex_mb2,-500E3,Model.rho_c,Model.rho_m,25E3,Model);
    V_flex(1,3) = 0; V_flex(2,2) = 0; V_flex(3,3) = 0;
    [~,DV_5] = degreeVariance(V_flex);
    rmse_array(i) = rmse(DV_0(2:90),DV_5(2:90));
end

figure(10)
plot(Te/1e3,rmse_array*10^10);
xlabel('Te (km)'); ylabel('RMS error (mGal^2)'); xlim([40 150]);
title("Root mean square error for various Te values")

[n,phi,phi2,flex_mb] = flexural_model(120e3,Model,airy_sc,lon,lat);

figure(11)
semilogx(n,phi)
hold on
semilogx(n,phi2)
xlim([2 180]); xticks([2 5 10 25 50 100])
xlabel('Spherical Harmonic Degree (n)'); ylabel('Flexural Response Function')
title('Flexural response function'); legend("Flat plate","Shell")
colormap(jet);
figure(12)
imagesc(lon,flip(lat),-flex_mb/10^3);c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)')
ylabel(c,'Crustal depth (km)'); title('Flexural model (shell)')
set(gca,'YDir','normal')
colormap(jet);

clear n phi phi2 Te temp temp_moho

%% Question 8 - Lateral Density Variations
% scaling down of gravity data
dataG = flip(grav_data.vec.R);
dataG = imresize(dataG,[45,90]);
dataG = imresize(dataG,[180,360]);

rho_c_min = 2400; % [kg/m^3]
Model.rho_c = 2700*ones(size(topo)); % [kg/m^3]
Model.moho = -flex_mb; % crust-mantle boundary

k = 0.2; % scaling factor

err = 1000;
iter = 1;
threshold = 200; % [mGal]
residual = 0;
while err>threshold && iter<31

    Model.rho_c = Model.rho_c+k*residual*10^5; % density adjustment
    % prevent crust density from being above mantle density
    for i = 1:size(Model.rho_c,1)
        for j = 1:size(Model.rho_c,2)
            if Model.rho_c(i,j) > Model.rho_m
                Model.rho_c(i,j) = Model.rho_m;
            end
        end
    end
    
    % V of model
    V_boug = segment_2layer_model(topo,Model.moho,-500E3,Model.rho_c,Model.rho_m,25E3,Model);
    V_boug(1,3) = 0; V_boug(2,2) = 0; V_boug(3,3) = 0;

    % gravity data of model
    [model_data] = model_SH_synthesis([0.5 359.5 1],[-89.5 89.5 1],0,[0 50],V_boug,Model);

    % residual and error calcuation
    residual = (dataG - flip(model_data.vec.R));
    if iter == 1
        residual0 = residual;
    end
    err = max(max(abs(residual)))*10^5;
    fprintf("Iteration %.0f error: %.2f; max density: %.0f.\n",iter,err,max(max(Model.rho_c)))

    iter = iter+1;
end

% Visualization
figure(13)
subplot(3,2,1) % observed gravity
imagesc(lon,flip(lat),dataG*10^5); c=colorbar; colormap(jet)
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('DataG'); set(gca,'YDir','normal')
subplot(3,2,2) % modeled gravity
imagesc(lon,flip(lat),flip(model_data.vec.R)*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('ModelG'); set(gca,'YDir','normal')
subplot(3,2,3) % topography
imagesc(lon,flip(lat),topo/1e3); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Topography (km)'); title('Topography'); set(gca,'YDir','normal')
subplot(3,2,4) % crustal density
imagesc(lon,flip(lat),Model.rho_c); c=colorbar; set(gca,'YDir','normal')
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Crustal density (kg/m^3)'); title('Crustal density');
c.Ticks = [2500 2750 3000 3250];
subplot(3,2,5)
imagesc(lon,flip(lat),residual0*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('Initial residual');
set(gca,'YDir','normal')
subplot(3,2,6)
imagesc(lon,flip(lat),residual*10^5); c=colorbar;
xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
ylabel(c,'Potential mGal'); title('Final residual')
set(gca,'YDir','normal');

% Lateral density inversion
figure(14)
axesm('MapProjection','mollweid', 'Frame', 'on', 'Grid', 'on');
temp = flipud(Model.rho_c); dens = zeros(size(Model.rho_c));
dens(:,1:180) = temp(:,181:360); dens(:,181:360) = temp(:,1:180);
geoshow(dens, [1 90 180], 'DisplayType', 'surface'); c=colorbar;
ylabel(c,'Crustal density (kg/m^3)');
colormap(jet);

clear i j k iter residual residual0 threshold err
%%
gmt = matrix2gmt(flex_mb,lon, lat);

%% Flexural Model Calculation
function[n,phi,phi2,flex_mb2] = flexural_model(Te,Model,airy_sc,lon,lat)
    E = 65*10^9; % [Pa] - source: Kalousova, Soucek, & Cadek, 2012
    v = 0.25; % [poisson's ratio] - source:Kalousova, Soucek, & Cadek, 2012
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

%% Vector converter
function vector=convert(A)
    v=zeros(180*360,1);
    i=1;
    for lats=1:1:180
        for lons=1:1:360
            v(i)=A(lats, lons);
            i=i+1;
        end
    end
    vector=v;
end