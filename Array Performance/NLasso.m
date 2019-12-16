function [omegaList, gainList, residueList] = NLasso(y,Phi,txf,tyf, ...
			      	   tau , R_c,omegaList,gainList)
% SUMMARY:
% 
%   given measurements: y = S * (mixture of sinusoids) + white noise
%          and a parameter tau which relates reduction in residue to
%          sparsity of the explanation (number of sinusoids), 
%          this module returns two **ordered** lists: a list of 
%          estimates of frequencies of the sinusoids in the mixture 
%          and and another list of corresponding gains 
% INPUT:
%    y - measurements
%    S - measurement matrix: can be compressive, or
%        just the identity matrix
%    tau - algorithm parameter which determines what the minimum payoff
%          we expect per sinusoid
%        - should be slightly more than the noise level sigma2
%        We solve for gains and continuous frequencies which minimize:
%        minimize norm(y - sum of sinusoids)^2 + tau * ell0_norm(gain)
%        ** LARGE VALUES OF TAU PROMOTE SPARSITY **% 
%    overSamplingRate (optional) -  used to determine how finely 
%              we sample frequency grid on which we precompute
%              responses when measurements are compressive (using 
%              IFFTs) how fine (in terms of multiples of the FFT 
%              grid) do we want the coarse grid to be?
%              number of grid points = overSamplingRate * N
%              Default value is 4
% 
%    R_s (optional) : number of Newton steps for a Single frequency update
%                           Default value is 1
%    R_c (optional) : number of rounds of Cyclic Refinement for all of the
%                           already detected sinusoids
%                           Default value is 3
% 
% OUTPUT:
%   omegaAns    - frequencies
%   gainAns     - gains of estimated frequencies
%   residueList - trajectory of the energy in the residual
%                 measurements as we add new sinsoids

VERBOSE = 0;

% if ~exist('overSamplingRate_1','var'), overSamplingRate_1 = 4;
% elseif isempty(overSamplingRate_1), overSamplingRate_1 = 4; end
% 
% if ~exist('Nref','var'), Nref = 4;
% elseif isempty(Nref), Nref = 4; end
R_s=2;
if ~exist('R_c','var'), R_c = 3;
elseif isempty(R_c),    R_c = 3; end

% Initial Setup
N1 = size(y , 1);
N2 = size(y , 2);

ant_idx_1 = txf - mean(txf); 
ant_idx_2 = tyf - mean(tyf); 

sinusoid_1d = @(w,ant_idx) exp(1j* ant_idx.'* w);
sinusoid_2d = @(w1,ant_idx_1,w2,ant_idx_2) Phi*(sinusoid_1d(w1,ant_idx_1).*sinusoid_1d(w2,ant_idx_2))/sqrt(length(ant_idx_1));

    [omegaList, gainList, y_r] = solveLeastSquares(y , Phi, omegaList, ...
         ant_idx_1, ant_idx_2); 

% Vectorize the signal
y_r_vec = y_r(:);

% Keeping track of residual energy
residueList = y_r_vec' * y_r_vec;

% while true
    
 
    % Add newly detected sinusoid to the ordered lists
%     omegaList = [omegaList; omega_new];
%     gainList  = [gainList; gain_new];
    
    % refine all frequencies detected so far
    % can be interpreted as a search for better frequency supports
    [omegaList, gainList, y_r] = refineAll(y_r, Phi, omegaList,...
        gainList, ant_idx_1, ant_idx_2 , R_s, R_c);
    % refineAll only uses refineOne to tweak parameters and the energy 
    % in the residual measurements y_r can only decrease as a result
    
    if (VERBOSE) fprintf('All sinusoids refined:\n'); omegaList 
    end;

    % Solve least squares for the dictionary set [Ax(omega)] omega in 
    % omegaList
    [omegaList, gainList, y_r] = solveLeastSquares(y , Phi, omegaList, ...
         ant_idx_1, ant_idx_2);    
    % ensures that for the support we have settled on our choice of 
    % gains is optimal (in terms giving us the lowest residual energy)
    % though the practical value of the least-squares step is debatable
    % when the frequency support is well-conditioned (since least-squares
    % is a by-product of refineAll once the frequencies have converged)
    % we need this step for theoretical guarantees on convergence rates
 
    residue_new = y_r(:)'*y_r(:);
    residueList = [residueList; residue_new];
        % Stopping rule:
%     if residue_new <= tau
%         break;
%     end   
% end

% revert to standard notion of sinusoid: 
%           exp(1j*(0:(N-1))'*omega)/sqrt(N)
% gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
% omegaList = wrap_2pi(omegaList);

end

% -----------------------------------------------------------------

function [omega, gain, y_r] = detectNew(y, Phi,...
					 ant_idx_1, ant_idx_2,posa)
% SUMMARY:
% 
% 	detects a new sinusoid on the coarse grid
% 
% INPUT:
% 	y - measurements
% 	sampledManifold - when measurements are compressive,
% 		  we precompute (using IFFT operation) a
% 		  **dictionary** of responses for sinusoids 
%  		  corresponding to frequencies on a coarse grid 
% 		  and store them in the MATLAB **structure** 
% 		  sampledManifold
%
% OUTPUT:
% 	omega - frequency on [0,2*pi) which best explains
% 			the measurements y
% 	gain  - corresponding complex gain
% 	y_r   - after removing the detected sinusoid from the
%         	measurements y, y_r is the ***residual measurement***
%          res_inf_normSq_rot - max energy among DFT directions - needed
%          for stopping criterion

% R = length(sampledManifold.coarseOmega);
% N = sampledManifold.length;
% OSR = round(R/N);

sinusoid_1d = @(w,ant_idx) exp(1j* ant_idx* w);
sinusoid_2d = @(w1,ant_idx_1,w2,ant_idx_2) Phi*(sinusoid_1d(w1,ant_idx_1).*sinusoid_1d(w2,ant_idx_2))/sqrt(length(ant_idx_1));

Dict=Phi*exp(1j*(2*pi*(ant_idx_1*(posa(1,:))+ant_idx_2*(posa(2,:)))));
% size(sum(abs(Dict).^2,1))
% size(sum(abs(Dict'*y).^2,2))
            proj_res=((sum(abs(Dict'*y).^2,2)).')./sum(abs(Dict).^2,1);
            [~,id1]=max(proj_res);
%             Dictf_shift=exp(1j*(2*pi*(SU(mask2(id1))*ant_idx_1+SV(mask2(id1))*ant_idx_2)));
%             Dictf_shifted=Phi*diag(Dictf_shift)*Dictf;
% proj_resf=sum(abs((Dictf_shifted')*y).^2,2);
% [~,id2]=max(proj_resf);

omega_1 = posa(1,id1)*2*pi;
omega_2 = posa(2,id1)*2*pi;
% 
% [~,id2]=max(proj_res);
% if 0
%     Dictf_shift=exp(1j*(2*pi*(posa3(1,id2)*ant_idx_1+posa3(2,id2)*ant_idx_2)));
%     Dictf_shifted=Phi*diag(Dictf_shift)*Dictf;
%     proj_resf=sum(abs((Dictf_shifted')*y).^2,2);
%     [~,id3]=max(proj_resf);
%     omega_1 = (posa3(1,id2)+SUf(id3))*2*pi;
%     omega_2 = (posa3(2,id2)+SVf(id3))*2*pi;
% else
%     omega_1 = (posa3(1,id2))*2*pi;
%     omega_2 = (posa3(2,id2))*2*pi;
% end
% compute the corresponding sinusoid:
x = sinusoid_2d(omega_1 , ant_idx_1, omega_2 , ant_idx_2);

% Note: gain = (x^{H}x)^{-1}*x^{H}*y; but x is unit norm, so:
gain = x(:)'*y(:); 
omega = [omega_1, omega_2];

% residual measurements after subtracting
% out the newly detected sinusoid
y_r = y - gain * x;
end

% --------------------------------------------------------------------

function omega_prime= wrap_2pi(omega)
% SUMMARY: Restricts frequencies to [0, 2*pi)
% INPUT: A vector of spatial frequencies
% OUTPUT: A vector of coresponding spatial frequencies in [0, 2*pi)

omega_prime = angle(exp(1j*omega));
omega_prime(omega_prime < 0) = omega_prime(omega_prime < 0) + 2*pi;

end

% --------------------------------------------------------------------

function [omega, gain, y_r] = refineOne(y_r,Phi, omega, gain, ...
			 ant_idx_1, ant_idx_2 , isOrth)
% SUMMARY:
%   Refines parameters (gain and frequency) of a single sinusoid
%   and updates the residual measurement vector y_r to reflect
%   the refinement -- This function applies one Newton step only.
% INPUT:
% 	y_r - residual measurement (all detected sinusoids removed)
%	omega - current estimate of frequency of sinusoid we want to
% 			refine
%	gain - current estimate of gain of sinusoid we want to refine
% 	S - measurement matrix
%           if [], measurements are direct i.e S = eye(N),
% 	    where N is length of sinusoid	
% 	ant_idx - translating indexes to phases in definition of sinusoid
%	isOrth - binary flag - is y_r orthogonal to x(omega) 
%	       - default - false
% OUTPUT:
%       refined versions of omega, gain and y_r
%       (see INPUT for definitions)

if ~exist('isOrth', 'var'), isOrth = false; end;

sinusoid_1d = @(w,ant_idx) exp(1j* ant_idx* w);
sinusoid_2d = @(w1,ant_idx_1,w2,ant_idx_2) (sinusoid_1d(w1,ant_idx_1).*sinusoid_1d(w2,ant_idx_2))/sqrt(length(ant_idx_1));

N1 = length(ant_idx_1);
N2 = length(ant_idx_2);
% 
% A = repmat(ant_idx_1(:), 1, N2);
% B = repmat(ant_idx_2(:).', N1, 1);

x_theta  = sinusoid_2d(omega(1) , ant_idx_1, omega(2) , ant_idx_2);

% dx_theta = 1j * ant_idx .* x_theta;
% d2x_theta = -(ant_idx.^2) .* x_theta;

dx_domega_1 = Phi*1j*(ant_idx_1.*x_theta);
dx_domega_2 = Phi*1j*(ant_idx_2.*x_theta);

d2x_domega_1 = -Phi*((ant_idx_1.^2).*x_theta);
d2x_domega_2 = -Phi*((ant_idx_2.^2).*x_theta);

% add the current estimate of the sinusoid to residue
y = y_r + gain*Phi*x_theta;


% UPDATE GAIN
% recompute gain and residue to ensure that 
% y_r is orthogonal to x_theta - this property
% is lost when we refine other sinusoids
if ~isOrth    
    gain = (Phi*x_theta(:))'*y(:);
    y_r = y - gain*Phi*x_theta;
end

% Compute derivatives
der1_omega_1 = -2*real(gain * y_r(:)'*dx_domega_1(:));
der2_omega_1 = -2*real(gain * y_r(:)'*d2x_domega_1(:)) +...
    2*abs(gain)^2*(dx_domega_1(:)'*dx_domega_1(:));

der1_omega_2 = -2*real(gain * y_r(:)'*dx_domega_2(:));
der2_omega_2 = -2*real(gain * y_r(:)'*d2x_domega_2(:)) +...
    2*abs(gain)^2*(dx_domega_2(:)'*dx_domega_2(:));

% UPDATE OMEGA(1)
if der2_omega_1 > 0
    omega_1_next = omega(1) - der1_omega_1/der2_omega_1;
else
    omega_1_next = omega(1) - sign(der1_omega_1)*(1/4)*(2*pi/N1)*rand(1);
end

% UPDATE OMEGA(2)
if der2_omega_2 > 0
    omega_2_next = omega(2) - der1_omega_2/der2_omega_2;
else
    omega_2_next = omega(2) - sign(der1_omega_2)*(1/4)*(2*pi/N2)*rand(1);
end


% COMPUTE x_theta for omega_next so that we can compute 
% gains_next and y_r_next
% x_theta  = exp(1j*ant_idx*omega_next)/sqrt(N);
x_theta = sinusoid_2d(omega_1_next , ant_idx_1, omega_2_next , ant_idx_2);


% UPDATE GAIN
gain_next = (Phi*x_theta(:))'*y(:);

% UPDATE RESIDUE
y_r_next = y - gain_next*Phi*x_theta;

% check for decrease in residue -  needed as a result of 
% non-convexity of residue (even when the cost surface 
% is locally convex); This is the same as checking whether 
% |<y, x_theta>|^2/<x_theta,x_theta> improves as we ensured 
% that y - gain*x_theta is perp to x_theta by recomputing gain

if (y_r_next(:)'*y_r_next(:)) <= (y_r(:)'*y_r(:))
    % commit only if the residue decreases
    omega = [omega_1_next, omega_2_next];
    gain = gain_next;
    y_r = y_r_next;
end

end

% --------------------------------------------------------------------

function [omegaList, gainList, y_r] = solveLeastSquares(y ,Phi, omegaList, ...
     ant_idx_1, ant_idx_2)
% SUMMARY:
%    Reestimates all gains wholesale
% N1 = length(ant_idx_1);
% N2 = length(ant_idx_2);

sinusoid_1d = @(w,ant_idx) exp(1j* ant_idx* w);
sinusoid_2d = @(w1,ant_idx_1,w2,ant_idx_2) Phi*(sinusoid_1d(w1,ant_idx_1).*sinusoid_1d(w2,ant_idx_2))/sqrt(length(ant_idx_1));

for k = 1:size(omegaList , 1)
    this_x = sinusoid_2d(omegaList(k,1), ant_idx_1 , omegaList(k,2), ant_idx_2);
    if k == 1
        A = this_x(:);
    else
        A = [A , this_x(:)];
    end    
end    

% update gains
gainList = (A'*A)\(A'*y(:));
% update residues
y_r_vec = y(:) - A*gainList;

y_r = y_r_vec;%reshape(y_r_vec, N1, N2);

% energy in the residual measurement y_r is guaranteed to not increase 
% as a result of this operation. Therefore, we refrain from checking 
% whether the energy in y_r has increased
end

% --------------------------------------------------------------------

function [omegaList, gainList, y_r] = refineAll(y_r, Phi, omegaList,...
    gainList, ant_idx_1, ant_idx_2, R_s, R_c)
% SUMMARY:
%   uses refineOne algorithm to refine frequencies and gains of
%   of all sinusoids
% INPUT:
%    y_r - residual measurement after all detected sinusoids have been
%          removed
%    omegaList - list of frequencies of known(detected) sinusoids
%    gainList  - list of gains of known(detected) sinusoids
%    S - measurement matrix
%        if [], measurements are direct i.e S = eye(N), where N is
%       length of sinusoid	
%    ant_idx - translating indexes to phases in definition of sinusoid
%    R_s - number of times each sinusoid in the list is refined
%    R_c - number of cycles of refinement fot all of frequencies that have
%       been estimated till now
%
% OUTPUT:
%       refined versions of inputs omegaList, gainList, y_r

K = size(omegaList , 1); % number of sinusoids


% Total rounds of cyclic refinement is "R_c"
for i = 1:R_c
    
    % chose an ordering for refinement
    
    % % *random* ordering
    % order = randperm(K);
    
    % *sequential* ordering
    order = 1:K;
    
    for j = 1:K
        l = order(j);
        
        
        
        % parameters of the l-th sinusoid
        omega = omegaList(l ,:);
        gain = gainList(l);
            
        % refinement repeated "R_s" times per sinusoid
        for kk = 1:R_s

            % refine our current estimates of (gain, omega) of the
            % l-th sinusoid
            [omega, gain, y_r] = refineOne(y_r,Phi,...
                       omega, gain, ant_idx_1, ant_idx_2, false);
                   
        end
        
        omegaList(l , :) = omega;
        gainList(l) = gain;
            % refineOne ensures that (gain, omega) pair are altered iff
            % the energy in the residual measurements y_r does not 
            % increase
        
    end
    
end

end

