function [ point , n_points,ninter] = get_intersection(X,Xs, u_n_intervs,v_n_intervs,u_knots,v_knots,...
     us_n_intervs,vs_n_intervs,us_knots,vs_knots,isec,n_isec, nit,nt )
 % computes intersection of BSpline patch with BSpline thread

n_points = 0;
ninter =  zeros(u_n_intervs,v_n_intervs);
for k=1:u_n_intervs
    us1 = u_knots(k+2);
    ue1 = u_knots(k+3);
    u1_c =(us1 + ue1)/2;
    ui = [us1 ue1] ;
    for l=1:v_n_intervs
        vs1 = v_knots(l+2);
        ve1 = v_knots(l+3);
        v1_c = (vs1  + ve1)/2;
        vi = [ vs1 ve1];
        s=0;
        if n_isec(k,l) ~= 0
            for p =1: n_isec(k,l)
                m = ceil(isec(k,l,p) / us_n_intervs);
                o = isec(k,l,p) - (m-1)*us_n_intervs;
                sp_i = (m-1)*(us_n_intervs) + o;
                % v2 fixed
                u2i = [us_knots(m+2) us_knots(m+3)];
                v_knot = linspace(vs_knots(o+2),vs_knots(o+3),nt);
                for h =1:length(v_knot)
                    u2_c = (us_knots(m+2) + us_knots(m+3))/2;
                    v2i = v_knot(h);
                    uvu2v2 = [u1_c;v1_c;u2_c]; % initial condition
                    [ uvu2v2,ptl,conv ] = patch_patch_intersection( uvu2v2, ui,vi,u2i,v2i, u_knots, v_knots, X,us_knots, vs_knots, Xs, nit ,m,o,k,l);
                    if conv ~= 0
                        s = s+1;
                        n_points = n_points +1;
                        point(n_points,:) = [uvu2v2',v_knot(h),ptl(1),ptl(2),ptl(3),k,l,m,o];
                    end
                end
                
                % u2 fixed
                v2i = [vs_knots(o+2) vs_knots(o+3)];
                u_knot = linspace(us_knots(m+2),us_knots(m+3),nt);
                for h =1:length(u_knot)
                    v2_c = (vs_knots(o+2) + vs_knots(o+3))/2;
                    u2i = u_knot(h);
                    uvu2v2 = [u1_c;v1_c;v2_c]; % initial condition
                    [ uvu2v2,ptl,conv ] = patch_patch_intersection( uvu2v2, ui,vi,u2i,v2i, u_knots, v_knots, X,us_knots, vs_knots, Xs, nit ,m,o,k,l);
                    if conv ~= 0
                        s = s+1;
                        n_points = n_points +1;
                        point(n_points,:) = [uvu2v2(1:2)',u_knot(h),uvu2v2(3),ptl(1),ptl(2),ptl(3),k,l,m,o];
                    end
                end
                ninter(k,l) = s;
            end
        end
    end
end


end



