function [a used_iter] = CoSaMP(z,Phi,s, max_iter)

% CoSaMP for sparse signal recovery
% inputs: s-->spasity
%         z-->measurements
%         Phi-->mensurement matrix
%         max_iter --> max iteration

% output: a-->the recovery signal

% CoSaMP: ACHA (2009) 301-321

d = size(Phi,2); % d-dimensinal signal


%Begin CoSaMP
err = 10^(-3);  % residual error
a = zeros(d,1); % Initialize the recovery signal
v = z;
it=0;
stop = 0;

while ~stop
    y = Phi'*v;
    [tmp, ix] = sort(abs(y), 'descend');
    Omega = ix(1:2*s);
    [tmp, ix] = sort(abs(a), 'descend');
    T = union(Omega, ix(1:s));
    %Signal Estimation
    b = zeros(d, 1);
    b(T) = Phi(:, T) \ z;
    %Prune
    [tmp, ix] = sort(abs(b), 'descend');
    a = zeros(d, 1);
    a(ix(1:s), 1) = b(ix(1:s), 1);
    %Sample Update
    v = z - Phi*a;
    %Iteration counter
    it = it + 1;
    
    %Check Halting Condition
    if (it > max_iter ||  norm(v)<=err*norm(z))   % norm(r_n)<=err  % Normally: max_iter = 6*(s+1)
      stop = 1;
    end
end %End CoSaMP iteration
used_iter = it;










