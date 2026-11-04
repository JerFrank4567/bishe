function [xr support support_num] = SAMP_distributed(z,Phi,step_size,TH,step_size_max,debug_en)
% SAMP: Sparsity Adaptive Matching Pursuit algoritm for compressed sensing.
% For theoretical analysis, please refer to the paper : 
% Thong. T. Do, Lu Gan and Trac D. Tran ,"Sparsity Adaptive Matching
% Purusit for practical compressed sensing" available at http://dsp.ece.rice.edu/cs
 
% Written by Thong Do(thongdo@jhu.edu)
% Updated on July, 26th 2008

% stop 条件改变为某一Block联合MMV分量的norm小于一定门限

% parameter usage:
%   z: Mx1 observation vector 
%   Phi: MxN measurement matrix
%   step_size: any positive integer value not larger than sparsity
%   sigma: noise energy when sensing
%   xr: reconstructed sparse signal
%   iter_num: number of iterations

m = size(Phi,2); % 
Block_len = m;
n = size(z,2); % m x n 的稀疏矩阵
a = zeros(m, n ); % 初始的稀疏矩阵

% Initialization
iter_num = 0; 
iter_num2 = 0; 
actset_size = step_size;
active_set = [];
res = z;
stg_idx = 1; % stage index
last_sparse_norm = +inf;
equ_err_flag=0;
while  1 

    % candidate list
    for pp = 1:n
        y(:,pp) = Phi(:,:,pp)'*res(:,pp);
    end
 
    [sort_index] = CoSaMP_MMV_Block_sub_sort(y,Block_len);
    candidate_set = union(active_set, sort_index(1:actset_size));
    
    b = zeros(m, n);
    [T_index] = CoSaMP_MMV_Block_sub_index(candidate_set,m,Block_len);
    for pp = 1:n
        b(T_index',pp) = Phi(:, T_index,pp) \ z(:,pp);
    end
    
    [sort_index2] = CoSaMP_MMV_Block_sub_sort(b,Block_len);
    active_set_new = sort_index2(1:actset_size);
    [T_active_set_new] = CoSaMP_MMV_Block_sub_index(active_set_new,m,Block_len);
    a = zeros(m, n ); 

    for pp = 1:n
        a(T_active_set_new,pp) = Phi(:, T_active_set_new,pp) \ z(:,pp);
    end
    
    %new_res = z - Phi*a;
    new_res = zeros(size(z));
    for pp = 1:n
        new_res(:,pp) = z(:,pp)-Phi(:,:,pp)*a(:,pp);
    end
    
%     % finalist
%     [val idx2] = sort(abs(pinv(Phi(:,candidate_set))*z), 'descend');
%     T_active_set_new = candidate_set(idx2(1:actset_size));
%     xr_active_set_temp = pinv(Phi(:,T_active_set_new))*z;
%      xr_est_temp(T_active_set_new)=xr_active_set_temp;
%     new_res = z-Phi(:,T_active_set_new)*xr_active_set_temp;
 
    % debug 
    a_temp = reshape(a,Block_len,[]);
    if debug_en==1
        res_norm(iter_num+1) = norm(new_res,'fro');

        stg_idx_record(iter_num+1) =  stg_idx;
        res_norm_sqr = res_norm.^2;
        
        figure(101);subplot(2,1,1);semilogy(res_norm_sqr);subplot(2,1,2);stem(stg_idx_record);
        figure(102);stem(sort_index2(1:actset_size), energy_avg( a_temp(sort_index2(1:actset_size),:)) );hold off;
    end

    if (norm(new_res,'fro') > last_sparse_norm)
%         iter_num2=0;
%         % shift into a new stage
%         %res_norm2(stg_idx)= norm(new_res);
%         
% %         last_active_set = active_set_new;
% %         T_active_set_last = T_active_set_new;
%         %break;
% %         last_sparse_norm = norm(new_res,'fro');
%             if stg_idx>1;
                
               T_active_set_last= T_active_set_last_stage;
              last_active_set= last_stage_active_set;
               break;
%             end
%             stg_idx = stg_idx + 1;
%         actset_size = stg_idx*step_size;
%         equ_err_flag=0;
    elseif (norm(new_res,'fro') >= norm(res,'fro'))
        iter_num2=0;
        % shift into a new stage
        %res_norm2(stg_idx)= norm(new_res);
        stg_idx = stg_idx + 1;
        actset_size = stg_idx*step_size;
        last_active_set = active_set_new;
        T_active_set_last = T_active_set_new;
        T_active_set_last_stage = T_active_set_new;
        last_stage_active_set = active_set_new;
        last_sparse_norm = norm(new_res,'fro');
        equ_err_flag=1;
    else
        
        % update residual and active set   
         iter_num2=iter_num2+1;
        res = new_res;
        last_active_set = active_set_new;
        active_set= active_set_new;
        T_active_set_last = T_active_set_new;
        equ_err_flag=0;
        
    end 
    
    
    iter_num = iter_num +1;
    
    
    temp_xxx =min(energy_avg( a_temp(active_set_new,:)));
     if  ( ( stg_idx>1 )&  temp_xxx<TH ) %|| (last_sparse_norm<norm(new_res,'fro')) %|| (stg_idx > step_size_max) || 
         %stg_idx>1 & iter_num2==1 & isempty( setdiff(last_active_set,active_set) ) & energy_avg( a_temp( setdiff(active_set,last_active_set),: ) )<TH % 最弱经能量低于门限，残差反增，超过一定稀疏度。
         T_active_set_last= T_active_set_last_stage;
        last_active_set= last_stage_active_set;
         break;
     end
    
    
end % looresp
 
% reconstruction
xr = zeros(m, n ); 
for pp=1:n
xr(T_active_set_last,pp) = pinv(Phi(:,T_active_set_last,pp))*z(:,pp);
end
support=sort(last_active_set,'ascend');
support_num = numel(support);