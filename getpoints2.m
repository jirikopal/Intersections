function [ X ] = getpoints2()
% function supply set of points X = [x_i,y_i,z_i,w_i] which defines surface
% x_i and y_i correspond to a equidistant grid  
% z_i is function value for corresponding [x_i,y_i], i.e., z_i = F(x_i,y_i)
% w_i is equal to 1
% X = [x,y,z,w]

nu = 20;
nv = 20;
h = 0;
X = zeros(nu*nv,4);

for i=1:nu
    for j = 1:nv
        h = h+1;
        x = 2*(i-1)/(nu-1);
        y = 2*(j-1)/(nv-1);
        X(h,1) = x;
        X(h,2) = y;
        %X(h,3) = cos(10*pi*x)*cos(3*pi*y)+exp(2*sin(x*y));
        %X(h,3) = cos(3*pi*exp(-x))*cos(3*pi*y^2)+exp(2*sin(x*y));
        %X(h,3) = 3*x*y;%-0.1; %- 0.5*rand(1,1) +10;%+exp(2*sin(x*y)); + 0.01*exp(x*y)+
        X(h,3) = sin(2 *x*y);
        X(h,4) = 1;
    end
end

end

