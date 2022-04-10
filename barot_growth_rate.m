function [growth] = barot_wn_midlat_p(lat,u)

% If you use this script in your research please cite our studies who used
% this script as well:

%Chemke R., Ming Y. & Yuval J. 2022. Nature Climate Change.


%This function calucaltes the growth rate as a fucntion of
%wavenumber (growth) by applying a linear normal mode instability analysis
%to the linearized absolute vorticity equation.

%input parameters:

%u - zonal wind
%lat - latitude


omega=7.292e-5;
radius=6371000;
y=(pi./180).*lat.*radius;
dy1=(y(2:end)+y(1:end-1))./2;
dy2=(dy1(2:end)-dy1(1:end-1));
dy=(y(2:end)-y(1:end-1)).*(cosd((lat(2:end)+lat(1:end-1))./2));



y1=(y(2)-y(1)).*(cosd((lat(2)+lat(1))./2));
y2=(y(end)-y(end-1)).*(cosd((lat(end)+lat(end-1))./2));

R=287.04;

    beta=2*omega*cosd(lat)./radius;

l1=1;
l2=length(lat);

   waven=0:1e-7:2e-5;
for k=1:length(waven)-1

for n=l1:l2
    
    if n==l2
        %mewaven(k)admim lo shel w
        A(n-l1+1,n-l1+1)=waven(k)*beta(n)-(waven(k)^2)*(waven(k)*u(n)) ...
        +(waven(k)*u(n))./(y2).*(-1./dy(n-1)-1./y2)...
        -(waven(k)./(y2)).*(((u(n-1)-u(n))./(dy(n-1)))-(u(n))./(y2));
        A(n-l1+1,n-l1+1-1)=((waven(k)*u(n))./(y2.*dy(n-1)));
        %mewaven(k)admim shel w
        B(n-l1+1,n-l1+1)=-(waven(k)^2)+(1./(y2)).*((-1./(dy(n-1)))-1./(y2));
        B(n-l1+1,n-l1+1-1)=(1)./(y2.*dy(n-1));

    elseif n==l1
        A(n-l1+1,n-l1+1)=waven(k)*beta(n)-(waven(k)^2)*(waven(k)*u(n)) ...
        +(waven(k)*u(n))./(y1).*(-1./y1-1./dy(n))...
        -(waven(k)./(y1)).*(((-u(n))./(y1))-(u(n)-u(n+1))./(dy(n)));
        A(n-l1+1,n-l1+1+1)=(waven(k)*u(n)./(y1.*dy(n)));
        B(n-l1+1,n-l1+1)=-(waven(k)^2)+(1./(y1)).*((-1./(y1))-1./(dy(n)));
        B(n-l1+1,n-l1+1+1)=(1)./(y1.*dy(n));
    elseif n>l1 && n<l2
    A(n-l1+1,n-l1+1-1)=(((waven(k)*u(n)))./(dy2(n-1).*dy(n-1)));
    A(n-l1+1,n-l1+1+1)=(((waven(k)*u(n)))./(dy2(n-1).*dy(n)));
    A(n-l1+1,n-l1+1)=waven(k)*beta(n)-(waven(k)^2)*(waven(k)*u(n)) ...
        +(waven(k)*u(n))./(dy2(n-1)).*(-1./dy(n-1)-1./dy(n))...
        -(waven(k)./(dy2(n-1))).*(((u(n-1)-u(n))./(dy(n-1)))-(u(n)-u(n+1))./(dy(n)));
    B(n-l1+1,n-l1+1-1)=(1)./(dy2(n-1).*dy(n-1));
    B(n-l1+1,n-l1+1+1)=(1)./(dy2(n-1).*dy(n));
    B(n-l1+1,n-l1+1)=-(waven(k)^2)+(1./(dy2(n-1))).*((-1./(dy(n-1)))-1./(dy(n)));
    %
    end
end
    A1=A(1:end,1:end);
    B1=B(1:end,1:end);
    if isempty(find(isnan(A1))) && isempty(find(isnan(B1)))
    wimag(:,k)=imag(eig(B1\A1));%AI*B
    wreal(:,k)=real(eig(B1\A1));
    else
    wimag(:,k)=nan;%AI*B
    wreal(:,k)=nan;
    end
%end
end
www=squeeze(max(wimag(:,:),[],1));


growth=www;

end