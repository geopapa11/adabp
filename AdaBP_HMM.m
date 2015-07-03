classdef AdaBP_HMM < handle
    % AdaBP_HMM is the class which creates the AdaBP structure and is
    % restricted to HMMs.
    % Author: geopapa
    % $ Date: 2014/01/22 18:12:21 $
    
   properties (GetAccess = private, SetAccess = private)
       T          % number of time points
       h          % information vector
       Jon        % Jon{k}       = J{k,k}      on-diagonal  terms
       Joff       % Joff{k}      = J{k,k+1}    off-diagonal terms
       msg_h_fwd  % msg_h_fwd{k} = h_{k->k+1}
       msg_J_fwd  % msg_J_fwd{k} = J_{k->k+1}
       msg_h_bwd  % msg_h_bwd{k} = h_{k->k-1}
       msg_J_bwd  % msg_J_bwd{k} = J_{k->k-1}
   end
   
   methods
       function obj = AdaBP_HMM(A, muV, Q, Y, mu1, cov1)
           obj.T         = size(Y,2);      obj.h         = cell(1,obj.T);
           obj.Jon       = cell(1,obj.T);  obj.Joff      = cell(1,obj.T);
           obj.msg_h_fwd = cell(1,obj.T);  obj.msg_J_fwd = cell(1,obj.T);
           obj.msg_h_bwd = cell(1,obj.T);  obj.msg_J_bwd = cell(1,obj.T);
           
           for k = 1:obj.T
               if k == 1
                   muV_kminus1 = mu1;
                   Q_kminus1   = cov1;
               else
                   if iscell(muV),           muV_kminus1 = muV{k-1};  else muV_kminus1 = muV;  end
                   if iscell(Q),             Q_kminus1   = Q{k-1};    else Q_kminus1   = Q;    end
                   if isempty(muV_kminus1),  muV_kminus1 = zeros(size(mu1));                   end
               end
               
               if k == obj.T
                   A_k   = zeros(size(Q_kminus1));
                   muV_k = zeros(size(muV_kminus1));
                   Q_k   = eye(size(Q_kminus1));
               else
                   if iscell(A),       A_k   = A{k};    else A_k   = A;    end
                   if iscell(muV),     muV_k = muV{k};  else muV_k = muV;  end
                   if iscell(Q),       Q_k   = Q{k};    else Q_k   = Q;    end
                   if isempty(muV_k),  muV_k = zeros(size(mu1));           end
               end
               
               tmp1        = A_k'/Q_k;
               tmp2        = tmp1*A_k;
               tmp2        = (tmp2 + tmp2')/2;
               obj.h{k}    = Q_kminus1\muV_kminus1 - tmp1*muV_k;
               obj.Jon{k}  = Q_kminus1\eye(size(Q_kminus1)) + tmp2;
               obj.Jon{k}  = (obj.Jon{k} + obj.Jon{k}')/2;
               obj.Joff{k} = -tmp1;
           end
           
           % Propagate forward
           propagate(obj,     1, obj.T);
                      
           % Propagate backward
           %propagate(obj, obj.T,     1);
           % <--- I am not using propagate function to evaluate the backward messages, 
           % because the standard implementation leads to numerical instabilities.
           % Note that all backward J messages must be NEGATIVE DEFINITE.
           for k = obj.T:-1:2
               if k == obj.T
                   msg_h_nghb = zeros(size(obj.h{k}));
                   msg_J_nghb = zeros(size(obj.Jon{k}));
               else
                   msg_h_nghb = obj.msg_h_bwd{k+1};
                   msg_J_nghb = obj.msg_J_bwd{k+1};
               end
               tmp1 = obj.Jon{k} + msg_J_nghb;  tmp1 = (tmp1 + tmp1')/2;
               tmp2 = -obj.Joff{k-1}/tmp1;
               obj.msg_h_bwd{k} = tmp2*(obj.h{k} + msg_h_nghb);
               
               if iscell(A),  A_kminus1 = A{k-1};  else  A_kminus1 = A;  end
               if iscell(Q),  Q_kminus1 = Q{k-1};  else  Q_kminus1 = Q;  end
               obj.msg_J_bwd{k} = -A_kminus1'/Q_kminus1*A_kminus1;
               obj.msg_J_bwd{k} = (obj.msg_J_bwd{k}+obj.msg_J_bwd{k}')/2;
           end
           
% % %            %% DO NOT ERASE THIS EVEN IF IT SEEMS OBSOLETE IN FUTURE VERSIONS.
% % %            %% IT IS VERY USEFUL FOR DEBUGGING PURPOSES, IN CASE EVER ANY OF 
% % %            %% THE J MESSAGES TURNS OUT TO BE POSITIVE DEFINITE
% % %            %% DEBUGGING MODE STARTS HERE
% % %            cnt1 = 1;
% % %            for k = 1:obj.T-1
% % %                ll = eig(obj.msg_J_fwd{k});
% % %                if any(ll>0) && any(ll<0)
% % %                    disp(['# ', num2str(k), ': forward']);
% % %                elseif all(ll<0)
% % %                    cnt1 = cnt1 + 1;
% % %                end
% % %            end
% % %            
% % %            cnt2 = 1;
% % %            for k = obj.T:-1:2
% % %                ll = eig(obj.msg_J_bwd{k});
% % %                if any(ll>0) && any(ll<0)
% % %                    disp(['# ', num2str(k), ': backward']);
% % %                    ll
% % %                elseif all(ll<0)
% % %                    cnt2 = cnt2 + 1;
% % %                end
% % %            end
% % %            disp([num2str(cnt1), '/', num2str(obj.T), ' forward messages are well-posed.']);
% % %            disp([num2str(cnt2), '/', num2str(obj.T), ' backward messages are well-posed.']);
% % %            %% DEBUGGING MODE FINISHES HERE
       end
       
       function update(obj, k, Y_k, C_k, muW_k, R_k)
           if isempty(muW_k),  muW_k = zeros(size(Y_k));  end
           tmp1       = C_k'/R_k;
           tmp2       = tmp1*C_k;
           tmp2       = (tmp2 + tmp2')/2;
           obj.h{k}   = obj.h{k}    + tmp1*(Y_k - muW_k);
           obj.Jon{k} = obj.Jon{k}  + tmp2;
           obj.Jon{k} = (obj.Jon{k} + obj.Jon{k}')/2;
       end
              
       function propagate(obj, wlk_cur, wlk_nxt)
           % Propagate forward
           if wlk_cur < wlk_nxt
               for k = wlk_cur:wlk_nxt-1
                   if k == 1
                       msg_h_nghb = zeros(size(obj.h{k}));
                       msg_J_nghb = zeros(size(obj.Jon{k}));
                   else
                       msg_h_nghb = obj.msg_h_fwd{k-1};
                       msg_J_nghb = obj.msg_J_fwd{k-1};
                   end
                   tmp1             = obj.Jon{k} + msg_J_nghb;
                   tmp1             = (tmp1 + tmp1')/2;
                   tmp2             = -obj.Joff{k}'/tmp1;
                   obj.msg_h_fwd{k} = tmp2*(obj.h{k} + msg_h_nghb);
                   obj.msg_J_fwd{k} = tmp2*obj.Joff{k};
                   obj.msg_J_fwd{k} = (obj.msg_J_fwd{k} + obj.msg_J_fwd{k}')/2;
               end
           % Propagate backward
           elseif wlk_cur > wlk_nxt
               for k = wlk_cur:-1:wlk_nxt+1
                   if k == obj.T
                       msg_h_nghb = zeros(size(obj.h{k}));
                       msg_J_nghb = zeros(size(obj.Jon{k}));
                   else
                       msg_h_nghb = obj.msg_h_bwd{k+1};
                       msg_J_nghb = obj.msg_J_bwd{k+1};
                   end
                   tmp1             = obj.Jon{k} + msg_J_nghb;
                   tmp1             = (tmp1 + tmp1')/2;
                   tmp2             = -obj.Joff{k-1}/tmp1;
                   obj.msg_h_bwd{k} = tmp2*(obj.h{k} + msg_h_nghb);
                   obj.msg_J_bwd{k} = tmp2*obj.Joff{k-1}';
                   obj.msg_J_bwd{k} = (obj.msg_J_bwd{k} + obj.msg_J_bwd{k}')/2;
               end
           end
       end
       
       function reset_msg(obj)
           obj.msg_h_fwd = cell(1,obj.T);  obj.msg_J_fwd = cell(1,obj.T);
           obj.msg_h_bwd = cell(1,obj.T);  obj.msg_J_bwd = cell(1,obj.T);
       end
       
       function [h, J] = eval_mrg(obj, k)
           if k == 1
               msg_h_left = zeros(size(obj.h{k}));  msg_J_left = zeros(size(obj.Jon{k}));
           else
               msg_h_left = obj.msg_h_fwd{k-1};     msg_J_left = obj.msg_J_fwd{k-1};
           end
           
           if k == obj.T
               msg_h_right = zeros(size(obj.h{k}));  msg_J_right = zeros(size(obj.Jon{k}));
           else
               msg_h_right = obj.msg_h_bwd{k+1};     msg_J_right = obj.msg_J_bwd{k+1};
           end
           
           h = obj.h{k}   + msg_h_left + msg_h_right;
           J = obj.Jon{k} + msg_J_left + msg_J_right;
       end
   end % methods
end % classdef