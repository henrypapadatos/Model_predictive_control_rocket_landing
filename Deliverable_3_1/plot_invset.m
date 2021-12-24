function [] = plot_invset(X_f,state_names,title,file_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   X_f           - polytope maximum invariant set 
%   state_names   - list of strings contraining the names of the state (e.g. ["vz", "z"])
%   title         - Suptitle of the graph
%   file_name     - name of the file under which the plot will be saved 
% OUTPUTS
%   Display a figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Plotting invariant sets")
% Get the dimension of the state 
[F,~] = double(X_f);
nx = size(F,2);

% Initialize the figure
fig = figure();
fig.Position = [100 100 1200 600];
sgtitle(title);

% Loop over the dimensions 
for dim1=1:nx-1
    subplot(1,nx-1,dim1)
    hold on
    % Plot the invariant set
    X_f.projection(dim1:dim1+1).plot("alpha",0.1,"color",'g');
    xlabel(state_names(dim1))
    ylabel(state_names(dim1+1))
end

% Eventually save the figure 
if nargin > 2
    saveas(fig,file_name);
end

end

