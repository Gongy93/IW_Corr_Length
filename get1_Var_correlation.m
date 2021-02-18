clc
clear
close all

%% Calculate the correlation length scale for internal waves 
% (see Mazloff et al. 2018, Gong et al. 2021)

% Aims:
% Based on the numerical modelling outputs of a representative variable
% (e.g. temperature, isotherm, velocity and etc.), we can characterize the
% spatial variability of internal waves.

% Input:
% Var [lon*lat*time]
% dx: horizontal resolution [unit: m]
% x0: selected box size, e.g. 20 km * 20 km
% x_interval = 5: do the loop every 5 grids (= 0.5 km)

% Output:
% Params of the correlation ellipse [lon*lat],
% including semimajor, semiminor axes and the ellipse orientation.

% Note:
% If the 3D variable matrix is big, this script might take long time.
% More efficient way is to set a larger value of x_interval (e.g. > 10),
% then the looping step would be quicker but the correlation ellipse
% has the lower resolution.

%  written by Gong Yankun in February 2021

%% Here I take the mode-2 amplitude as an example
load('Amp2_time.mat','Amp2_time');
Var = Amp2_time; clear Amp2_time % Dimension: [lon*lat*time]

% Params setup
dx = 100; % horizontal resolution [unit: m]
x0 = 20000; % selected box size: 20 km * 20 km

% Define the starting point and ending point to do the loop
start_lon = x0/dx; 
end_lon = size(Var,1) - x0/dx;
start_lat = x0/dx;
end_lat = size(Var,2) - x0/dx;

%% Do the loop to calculate the correlation length at each point (Var1)
x_interval = 5; % do the loop every 5 grids (= 0.5 km)

for lo = start_lon :x_interval: end_lon
   for la = start_lat :x_interval: end_lat
   
 Var1 = squeeze(Var(lo,la,:));     
% Reshape the variable in the small box
 Var2 = reshape(Var(lo-(x0/dx-x_interval):lo+(x0/dx-x_interval),...
       la-(x0/dx-x_interval):la+(x0/dx-x_interval),:),...
       [(x0/dx*2-x_interval*2+1)*(x0/dx*2-x_interval*2+1) size(Var,3)])';
   
% Calculate the correlation coefficient between Var1 and Var2 
 [rho,pval]= corr(Var1,Var2); 
 corr_Var = reshape(rho,[(x0/dx*2-10+1) (x0/dx*2-10+1)]); clear rho pval
 delta_lon = -(x0/dx-x_interval):(x0/dx-x_interval); 
 delta_lat = -(x0/dx-x_interval):(x0/dx-x_interval);

 %% Define the correlation bands
 corr_min = 0.75;
 corr_max = 0.85; 
 
 [a,b] = find(corr_Var>corr_min & corr_Var<corr_max);
   
 delta_x = delta_lon(a);
 delta_y = delta_lat(b);
   
 for m  = 1:length(a)
       z(m) = log(corr_Var(a(m),b(m)));
 end

 %% Calculate the parameters for correlation ellipses 
 a1 = zeros(size(Var,1),size(Var,2));
 a2 = zeros(size(Var,1),size(Var,2));
 a3 = zeros(size(Var,1),size(Var,2));
 
 % if the samples are less than 3, it means the variable can fit to a
 % ellipse in the box. We thus state the value NAN. 
  if length(a) < 3   
     a1(lo,la) = nan; a2(lo,la) = nan; a3(lo,la) = nan;      
  else
 
 % fit the parameter to a ellipse function in the small box (x0*x0)
  param = z / [delta_x.^2;delta_x.*delta_y;delta_y.^2];
  a1(lo,la) = param(1);
  a2(lo,la) = param(2);
  a3(lo,la) = param(3);
  
  end 
   
  clear delta_x delta_y z param a b param

   end
  disp(['--------------  Have completed ' ...
      num2str((lo-start_lon)/end_lon*100,'%3.1f') ' %  ---------------']);
  disp(datestr(clock))
end

%% semimajor, semiminor axes and the ellipse orientation
sida = 0.5*atand(a2./(a1-a3));
max_a = abs(sqrt(1./(a1.*cosd(sida).^2 + a2.*sind(sida).*cosd(sida) + a3.*sind(sida).^2)));
min_a = abs(sqrt(1./(a1.*sind(sida).^2 - a2.*sind(sida).*cosd(sida) + a3.*cosd(sida).^2)));

% Save the correlation ellipse parameter for plotting at the next step
save('Corr_ellipse_params.mat','-v7.3','a1','a2','a3','sida','max_a','min_a');


