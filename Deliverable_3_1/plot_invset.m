function [] = plot_invset(X_f,title_graph,X_lqr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   X_f           - polytope maximum invariant set 
%   title         - title of the graph
%   X_lqr         - Initial constrained set (here with LQR controller)
% OUTPUTS
%   Display a figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Plotting invariant sets")
[F,f] = double(X_f);
nx = size(F,2);
for dim1=1:nx-1
    subplot(1,nx-1,dim1)
    hold on
    % Question to ask to TA
    if nargin > 2
        X_lqr.projection(dim1:dim1+1).plot("alpha",0.2,"color",'r');
    end
    X_f.projection(dim1:dim1+1).plot("alpha",0.1,"color",'g');
    xlabel("dimension " + int2str(dim1))
    ylabel("dimension " + int2str(dim1+1))
end
end

