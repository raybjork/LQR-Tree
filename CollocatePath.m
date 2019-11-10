function [trajectory] = CollocatePath(start, finish, dynamics, n)
    % start and finish should be state space positions
    % n is the number of knot points in collocation
    % dynamics should be function handle that takes state vector, and
    % return xdot
    %trajectory will be in form u1,x2,u2,dt and assume
    %linear splines for each u between knot points.
    
    %need to index decision variables, then construct dynamic constraints (
    
    
end

