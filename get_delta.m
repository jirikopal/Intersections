function [ dist ] = get_delta(u_knots, v_knots,us_knots, vs_knots,uvt,u2i,v2i,k,l,m,o,X,Xs)
% returns distance in x,y,z space 

  [uf, ~] = splinebasevec(u_knots,uvt(1),0,k);
    [vf, ~] = splinebasevec(v_knots,uvt(2),0,l);
   
    if length(u2i) == 1
        [v2f, ~] = splinebasevec(vs_knots,uvt(3),0,o);
        [u2f, ~] = splinebasevec(us_knots,u2i,0,m);
    end
    if length(v2i) == 1
        [u2f, ~] = splinebasevec(us_knots,uvt(3),0,m);
        [v2f, ~] = splinebasevec(vs_knots,v2i,0,o);
    end  
    
    deltaXYZ = (kron(vf',uf') * X)' - (kron(v2f',u2f') * Xs)';
    dist = norm(deltaXYZ);

end

