function [max_growth_midlat, most_us_wn, growth] = baroc_growth_rate(u,lon,theta,z,lat, up_b,t_mean)

% If you use this script in your research please cite one/more of our studies who used
% this script as well:

%growth rate papers:
%Chemke R., Ming Y. & Yuval J. (2022)
%Chemke R. (2022) Nature Communications
%Chemke R., Zanna L., Polvani L.M., Orbe C. & Sentman L. (2022) Journal of Climate
%Chemke R. & Ming Y. (2020) Geophysical Research Letters
%Chemke R. & Polvani L. M. (2019) Journal of Climate

%most unstable wn papers:
%Chemke R. (2017) Quarterly Journal of the Royal Meteorological Society
%Chemke R., Dror T. & Kaspi Y. (2016) Geophysical Research Letters
%Chemke R. & Kaspi Y. (2016) Journal Of The Atmospheric Sciences
%Chemke R. & Kaspi Y. (2016) Geophysical Research Letters.
%Chemke R. & Kaspi Y. (2015) Journal Of The Atmospheric Sciences
%


%This function calucaltes the maximum growth rate (max_growth_midlat), most
%unstable wavenumber (most_us_wn) and the growth rate as a fucntion of
%wavenumber (growth) by applying a linear normal mode instability analysis
%to the linearized QG equations.

%input parameters:

%u - zonal wind
%lon - longitude
%theta - potential temperature
%z - pressure
%lat - latitude
%up_b - upper boundary (tropopuase height)
%t_mean - temperature

%here the meridional wavenumber (l) is zero

omega=7.292e-5;
radius=6371000;
h1=up_b;


dz1=(z(2:end)+z(1:end-1))./2;%half levels
dz2(2:length(dz1))=(dz1(2:end)-dz1(1:end-1));
dz2(1)=(z(2)-z(1));
dz2(end+1)=(z(end)-z(end-1));
w=0:round((length(lon)-1)/3);
R=287.04;

theta_mean=theta;
h2=length(z);
f=2*omega*sind(lat);
    beta=2*omega*cosd(lat)./radius;

rho=z./(R.*t_mean);
Nz=-(1./(rho(1:end-1).*theta_mean(1:end-1))).*(theta_mean(2:end)-theta_mean(1:end-1));


   % waven=linspace(5*10^-7,5*10^-6,200);% for local NA
    waven=2*pi*w./(abs(lon(1)-lon(end))/180*pi*radius*cosd(lat));%for zonal mean low res.
    waven=0:1e-7:2e-5;%for zonal mean high res.
    waven1=waven;
    waven1(:)=0;
    v=zeros(size(u,1));
for k=1:length(waven)-1
%for l=1%;%1:length(waven1) %
l=1;
for n=h1:h2
    
    if n==h2
        A(n-h1+1,n-h1+1)=waven(k)*beta-(waven(k)^2+waven1(l)^2)*(waven(k)*u(n)+waven1(l)*v(n)) ...
               +((waven(k)*u(n)+waven1(l)*v(n))*(f.^2)./(dz2(n))).*(-1./(Nz(n-1)))...
              -(waven(k)*(f.^2)./(dz2(n))).*((u(n-1)-u(n))./(Nz(n-1)))...
              +(waven1(l)*(f.^2)./(dz2(n))).*((v(n-1)-v(n))./(Nz(n-1)));
        A(n-h1+1,n-h1+1-1)=((waven(k)*u(n)+waven1(l)*v(n))*(f.^2)./(dz2(n).*Nz(n-1)));
        B(n-h1+1,n-h1+1)=-(waven(k)^2+waven1(l)^2)+((f.^2)./(dz2(n))).*(-1./(Nz(n-1)));
        B(n-h1+1,n-h1+1-1)=(f.^2)./(dz2(n).*Nz(n-1));

    elseif n==h1
        A(n-h1+1,n-h1+1)=waven(k)*beta-(waven(k)^2+waven1(l)^2)*(waven(k)*u(n)+waven1(l)*v(n)) ...
            +((waven(k)*u(n)+waven1(l)*v(n))*(f.^2)./(dz2(n))).*((-1./(Nz(n))))...
            -(waven(k)*(f.^2)./(dz2(n))).*(((u(n+1)-u(n))./(Nz(n))))...
            +(waven1(l)*(f.^2)./(dz2(n))).*(((v(n+1)-v(n))./(Nz(n))));
        A(n-h1+1,n-h1+1+1)=(waven(k)*u(n)*(f.^2)./(dz2(n).*Nz(n)));
        B(n-h1+1,n-h1+1)=-(waven(k)^2+waven1(l)^2)+((f.^2)./(dz2(n))).*((-1./(Nz(n))));
        B(n-h1+1,n-h1+1+1)=(f.^2)./(dz2(n).*Nz(n));
    elseif n>h1 && n<h2
    A(n-h1+1,n-h1+1-1)=(((waven(k)*u(n)+waven1(l)*v(n)))*(f.^2)./(dz2(n).*Nz(n-1)));
    A(n-h1+1,n-h1+1+1)=(((waven(k)*u(n)+waven1(l)*v(n)))*(f.^2)./(dz2(n).*Nz(n)));
    A(n-h1+1,n-h1+1)=waven(k)*beta-(waven(k)^2+waven1(l)^2)*(waven(k)*u(n)+waven1(l)*v(n)) ...
        +((waven(k)*u(n)+waven1(l)*v(n))*(f.^2)./(dz2(n))).*((-1./(Nz(n-1)))-1./(Nz(n)))...
        -(waven(k)*(f.^2)./(dz2(n))).*(((u(n-1)-u(n))./(Nz(n-1)))-(u(n)-u(n+1))./(Nz(n)))...
        +(waven1(l)*(f.^2)./(dz2(n))).*(((v(n-1)-v(n))./(Nz(n-1)))-(v(n)-v(n+1))./(Nz(n)));
    B(n-h1+1,n-h1+1-1)=(f.^2)./(dz2(n).*Nz(n-1));
    B(n-h1+1,n-h1+1+1)=(f.^2)./(dz2(n).*Nz(n));
    B(n-h1+1,n-h1+1)=-(waven(k)^2+waven1(l)^2)+((f.^2)./(dz2(n))).*((-1./(Nz(n-1)))-1./(Nz(n)));
    
    end
end
    A1=A(1:end,1:end);
    B1=B(1:end,1:end);
    if isempty(find(isnan(A1))) && isempty(find(isnan(B1)))
    wimag(:,k)=imag(eig(B1\A1));
    wreal(:,k)=real(eig(B1\A1));
    else
    wimag(:,k)=nan;
    wreal(:,k)=nan;
    end
%end
end
www=squeeze(max(wimag(:,:),[],1));

gg=find(www(1:end-1)-www(2:end)>0);
if isempty(gg)==0
k_max=find(max(www(1:gg(end)+1))==www(1:gg(end)+1));
baroc_wn_midlat1=max(www(1:gg(end)+1));
else
k_max=nan;
baroc_wn_midlat1=nan;
end


most_us_wn=2*pi*radius*cosd(lat)./k_max;
max_growth_midlat=baroc_wn_midlat1;
growth=www;

end