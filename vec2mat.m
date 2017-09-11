function [ Z_coor ] = vec2mat( zs,u_n_basf,v_n_basf )
%compose matrix from vectors
Z_coor = zeros(u_n_basf,v_n_basf);

for i =1:v_n_basf
    Z_coor(:,i) = zs(1+(i-1)*(u_n_basf):i*(u_n_basf));
end

end

