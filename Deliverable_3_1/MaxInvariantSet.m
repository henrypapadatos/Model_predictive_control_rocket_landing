function [O_inf] = MaxInvariantSet(X,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   X           - polyhedron defined by X and U constraints
%   A           - matrix A such that x+ = Ax (for LQR A := (A + B*K))
% OUTPUTS
%   O_inf       - polytope Maximum invariant set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the maximal invariant set
i = 1;
O_inf = X;
while 1
	Oprev = O_inf;
	[F,f] = double(O_inf);	
    disp(size(F))
	% Compute the pre-set
	O_inf = polytope([F;F*A],[f;f]);
    if O_inf == Oprev 
        break; 
    end
	fprintf('Iteration %i... not yet equal\n', i)
	i = i + 1;
end

end

