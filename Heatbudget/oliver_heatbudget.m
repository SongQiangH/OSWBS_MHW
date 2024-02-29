clc;clear;

xw = 140; xe = 175; ys = 50; yn = 60; 
h_data = ncread('/public/home/songqh/my_data/ECCO2/MXLDEPTH/MXLDEPTH.2022.nc','MXLDEPTH');
h = squeeze(nanmean(nanmean(h_data(540:720,540:620,:),1),2));

R = 6371*1000; %m, earth radius specified in MOM4p1 (noted in mom4p0_manual)
A = R*R*(xe - xw)*(2*pi/360)*(sind(yn) - sind(ys)); %m^2, area of domain
load /public/home/songqh/my_data/project/MHW_2022/ECCO2_grid_bounds.mat
load /public/home/songqh/my_data/project/MHW_2022/ECCO2_grid_data_NWP.mat
DZ = diff(st_bnds_u);
nz = length(st_bnds_u);
xu = ECCO2_GRID.xu; yu = ECCO2_GRID.yu;
xt = ECCO2_GRID.xt; yt = ECCO2_GRID.yt; 
iyn = find(yu(:,1)-0.125==yn); iys = find(yu(:,1)-0.125==ys);
ixe = find(xu(1,:)-0.125 == xe); ixw = find(xu(1,:)-0.125 == xw); 
if isnan(iyn) == 1; display('yn not in u-grid'); end; 
if isnan(iys) == 1; display('ys not in u-grid'); end; 
if isnan(ixw) == 1; display('xw not in u-grid'); end; 
if isnan(ixe) == 1; display('xe not in u-grid'); end; 
IXt = ixw+1:ixe;  IYt = iys+1:iyn; %new grid
XT = xt(IYt,IXt); YT = yt(IYt,IXt); 
[ny,nx] = size(XT); 
for iy = 1:ny; for ix = 1:nx;
  dA(iy,ix) = R*R*0.25*(2*pi/360)*(sind(YT(iy,ix) + 0.125) - sind(YT(iy,ix) - 0.125));  
end; end;
dy_w = R*(0.25*2*pi/360); %m, meridional length
dy_e = R*(0.25*2*pi/360); %m, meridional length
dx_n = R*cosd(yn)*(0.25*2*pi/360); %m, zonal width 
dx_s = R*cosd(ys)*(0.25*2*pi/360); %m, zonal width 
int_dy_w = nansum(ny*dy_w); int_dy_e = nansum(ny*dy_e); 
int_dx_n = nansum(nx*dx_n); int_dx_s = nansum(nx*dx_s); 
%%
ncfile_sst = '/public/home/songqh/my_data/SST/oisst.2022_box.nc';
ncfile_u = '/public/home/songqh/my_data/ECCO2/UVEL/UVEL_interp/uvel.2022_box.nc';  
ncfile_v = '/public/home/songqh/my_data/ECCO2/VVEL/VVEL_interp/vvel.2022_box.nc';
ncfile_w = '/public/home/songqh/my_data/ECCO2/WVEL/WVEL_interp/wvel.2022_box.nc';
ncfile_theta = '/public/home/songqh/my_data/ECCO2/Theta/THETA_interp/THETA.2022_box.nc';
ncfile_q = '/public/home/songqh/my_data/ECCO2/Qnet/oceQnet.2022_box.nc';
ncfile_sw = '/public/home/songqh/my_data/ECCO2/Qsw/oceQsw.2022_box.nc';
u_data = ncread(ncfile_u, 'UVEL');
v_data = ncread(ncfile_v,'VVEL');
w_data = ncread(ncfile_w,'WVEL');
theta_data = ncread(ncfile_theta,'THETA');
sst_data = ncread(ncfile_sst, 'sst');
q_data = ncread(ncfile_q,'oceQnet');
sw_data = ncread(ncfile_sw,'oceQsw');
time = ncread(ncfile_u,'TIME');
nt = length(time); 
%%
temp = sst_data;
ue1 = squeeze(u_data(ixe,iys:iyn,1:2,:));
uw1 = squeeze(u_data(ixw,iys:iyn,1:2,:));
vn1 = squeeze(v_data(ixw:ixe,iyn,1:2,:));
vs1 = squeeze(v_data(ixw:ixe,iys,1:2,:));
q_data1 = squeeze(q_data(ixw:ixe-1,iys:iyn-1,:));
sw_data1 = squeeze(sw_data(ixw:ixe-1,iys:iyn-1,:));

temp = permute(temp, [3,2,1]);
ue1 = permute(ue1,[3,2,1]);
uw1 = permute(uw1,[3,2,1]);
vn1 = permute(vn1,[3,2,1]);
vs1 = permute(vs1,[3,2,1]);
q_data1 = permute(q_data1,[3,2,1]);
sw_data1 = permute(sw_data1,[3,2,1]);
%%
for iy = 1:length(iys:iyn - 1) %number of T grid cells in the meridional direction
  ue(:,:,iy) = nanmean(ue1(:,:,iy:iy+1),3); 
  uw(:,:,iy) = nanmean(uw1(:,:,iy:iy+1),3); 
end; 
for ix = 1:length(ixw:ixe - 1) %number of T grid cells in the zonal direction
  vn(:,:,ix) = nanmean(vn1(:,:,ix:ix+1),3); 
  vs(:,:,ix) = nanmean(vs1(:,:,ix:ix+1),3); 
end; 
tw = squeeze(nanmean(temp(:,iys+1:iyn,ixw:ixw+1),3)); 
te = squeeze(nanmean(temp(:,iys+1:iyn,ixe:ixe+1),3)); 
ts = squeeze(nanmean(temp(:,iys:iys+1,ixw+1:ixe),2)); 
tn = squeeze(nanmean(temp(:,iyn:iyn+1,ixw+1:ixe),2)); 
%%
tw(isnan(tw)) = 0;
te(isnan(te)) = 0;
ts(isnan(ts)) = 0;
tn(isnan(tn)) = 0;
%%
for it = 1:nt
  utw = squeeze(nanmean(uw(it,:,:),2))'.*squeeze(tw(it,:)); 
  int_utw(it) = nansum(utw*dy_w,2); %deg C x m^2/s; nt x (depth over top 60 m) 
  %
  ute = squeeze(nanmean(ue(it,:,:),2))'.*squeeze(te(it,:)); 
  int_ute(it) = nansum(ute*dy_e,2); %nt x (depth over top 60 m) 
%   %
  vtn = squeeze(nanmean(vn(it,:,:),2))'.*squeeze(tn(it,:));  
  int_vtn(it) = nansum(vtn*dx_n,2); %nt x (depth over top 60 m) 
%   %
  vts = squeeze(nanmean(vs(it,:,:),2))'.*squeeze(ts(it,:));  
  int_vts(it) = nansum(vts*dx_s,2); %nt x (depth over top 60 m) 
end; 
%%
int_utw_d = sum(int_utw'.*(ones(nt,2).*DZ),2)/sum(DZ); %deg C x m^2/s; nt (number of days)
int_ute_d = sum(int_ute'.*(ones(nt,2).*DZ),2)/sum(DZ); 
int_vtn_d = sum(int_vtn'.*(ones(nt,2).*DZ),2)/sum(DZ); 
int_vts_d = sum(int_vts'.*(ones(nt,2).*DZ),2)/sum(DZ); 
int_adv_d = (-(1/A)*(int_ute_d - int_utw_d + int_vtn_d - int_vts_d))*24*3600;
int_adv_d(isnan(int_adv_d)) = 0;
%%

% end; %it
cp = 3990; %J x deg C^-1 x kg^-1, specific heat at constant pressure
rho_ref = 1035; %kg/m^3, reference density applied in MOM4p1
for it = 1:nt; 
Qnet(it) = nanmean(nanmean(q_data1(it,:,:),2),3)/(h(it)*cp*rho_ref)*24*3600;
Qsw(it) = nanmean(nanmean(sw_data1(it,:,:),2),3)/(h(it)*cp*rho_ref)*24*3600; %deg C x s^-1 
end;
dTdt = diff(squeeze(nanmean(nanmean(temp(:,iys:iyn-1,ixw:ixe-1),2),3)));
%%
filename1 = '/public/home/songqh/my_data/project/MHW_2022/advect.mat';
filename2 = '/public/home/songqh/my_data/project/MHW_2022/Qnet.mat';
filename3 = '/public/home/songqh/my_data/project/MHW_2022/Qsw.mat';
filename4 = '/public/home/songqh/my_data/project/MHW_2022/dTdt.mat';

save(filename1, 'int_adv_d');
save(filename2, 'Qnet');
save(filename3, 'Qsw');
save(filename4, 'dTdt');