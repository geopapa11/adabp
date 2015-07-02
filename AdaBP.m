classdef AdaBP < handle
    % ADABP is the class which creates the AdaBP structure and returns
    % the (directed) message schedule between two nodes in time linear to
    % the length of the path between these two nodes.
    % Author: geopapa
    % $ Date: 2014/06/26 19:51:34 $

   properties (GetAccess = private, SetAccess = private)
       Ncc            % number of connected components
       N              % number of nodes
       Ne             % number of edges
       card           % cardinality (alphabet) of each node
       E              % edges
       pa             % parents
       ch             % children
       nghb           % neighbors (nodes take the original (absolute) indices)
       edgeOrig       % array of edge origins (left nodes of a node pair representing an edge)
       edge2idx       % map of edges to linear indices in edge matrix E
       nghb2idx       % map of neighbors to linear indices in cell array nghb
       sameCCNodes    % each cell array element contains the indices of nodes contained in that connected component
       node2compNode  % map of node's index (in the whole graph) to its index in the component it belongs to
       node2compID    % map of node's index (in the whole graph) to the component it belongs to
       prev_w         % component index of the node with the most recent measurement (from this component) 
       root           % root
       Eul            % euler tour
       De             % depth of each node in the Euler tour
       H              % index of the first occurence of node i in Eul
       I              % stores the index of the minimum element in a specified range
       isGaussian     % check if the graph is a Gaussian MRF
       isSumProd      % check if sum-product (goal: node marginals) or max-product (goal: map sequence) is used
       queue          % queue of messages to be updated at each iteration
       w_cur          % current walk element
       seq_prv        % MAP sequence at previous iteration
       lbl2idx        % Node labels to absolute indices
   end

   properties (GetAccess = private, SetAccess = private)
       phi      % node potentials
       psi      % edge potentials
       msg      % messages in Discrete Belief Propagation
       msg_idx  % indicator messages in Discrete Belief Propagation (used only in Max-Product form)
   end
   
   properties (GetAccess = private, SetAccess = private)
       h      % full information vector
       J      % full precision matrix
       msg_h  % h messages in Gaussian Belief Propagation
       msg_J  % J messages in Gaussian Belief Propagation
   end
   
   methods (Access = private)
       % Initialize graph; neighbors, children, parents, euler tour etc.
       function init(obj)
           % Propagate messages across all directions and determine neighbors
           if obj.isSumProd
               if obj.isGaussian
                   [~, ~, ~, ~, obj.msg_h, obj.msg_J, obj.nghb] = bp_ga(obj.h, obj.J);
               else
                   [~, obj.msg, obj.nghb] = bp(obj.E, obj.phi, obj.psi);
               end
           else
               [obj.seq_prv, ~, ~, obj.msg, obj.msg_idx, obj.nghb] = bp_maxsum(obj.E, obj.phi, obj.psi);
               obj.card  = cellfun(@length,obj.phi);
               obj.w_cur = 0;
               obj.queue = [];
           end
           
           % Get the number of edges in the graph
           obj.Ne = size(obj.E,1);
           
           % Create a map of edges to their linear indices in matrix E
           keys    = unique(obj.E(:,1));
           edgeIdx = 1:size(obj.E,1);
           obj.edgeOrig = sparse(keys,1,1:length(keys),max(obj.E(:)),1);
           obj.edge2idx = cell(length(keys),1);
           for k = 1:length(obj.edge2idx)
               idx = obj.E(:,1)==keys(k);
               obj.edge2idx{k} = sparse(obj.E(idx,2),1,edgeIdx(idx),max(obj.E(idx,2)),1);
           end
           
           % Create a map of neighbors to their linear indices in cell matrix nghb
           obj.nghb2idx = cell(obj.N,1);
           for i = 1:obj.N
               keys = obj.nghb{i};
               vals = 1:length(obj.nghb{i});
               if isempty(keys),  keys = 0;  vals = 0;  end
               obj.nghb2idx{i} = sparse(keys,1,vals,max(keys),1);
           end
           
           % Determine the connected components of the graph
           comp              = GraphUtils.conn_comp(obj.nghb);
           obj.Ncc           = length(comp);
           obj.sameCCNodes   = comp;
           obj.node2compNode = zeros(1,obj.N);
           obj.node2compID   = zeros(1,obj.N);
           comp_nghb         = cell(1,obj.Ncc);
           for i = 1:obj.Ncc
               obj.node2compNode(comp{i}) = 1:length(comp{i});
               obj.node2compID(comp{i})   = i;
               comp_nghb{i}               = obj.nghb(comp{i});
               if obj.Ncc > 1
                   for j = 1:length(comp_nghb{i})
                       for k = 1:length(comp_nghb{i}{j})
                           comp_nghb{i}{j}(k) = obj.node2compNode(comp_nghb{i}{j}(k));
                       end
                   end
               end
           end
           
           % Get the nodes' labels for each component
           obj.lbl2idx = cell(1,obj.Ncc);
           for i = 1:obj.Ncc
               lbl            = unique([comp_nghb{i}{:}]);
               obj.lbl2idx{i} = sparse(lbl,1,1:length(comp_nghb{i}),max(lbl),1);
           end
           
           obj.prev_w = zeros(1,obj.Ncc);
           
           obj.ch  = cell(1,obj.Ncc);  obj.pa = cell(1,obj.Ncc);
           obj.Eul = cell(1,obj.Ncc);  obj.De = cell(1,obj.Ncc);
           obj.H   = cell(1,obj.Ncc);  obj.I  = cell(1,obj.Ncc);
           for i = 1:obj.Ncc
               root_comp = randi(length(comp_nghb{i}));
               
               % Determine the children
               if obj.Ncc==1
                   obj.ch{1} = GraphUtils.get_ch(comp_nghb{1}, obj.root);
               else
                   obj.ch{i} = GraphUtils.get_ch(comp_nghb{i}, root_comp);
               end
               
               % Determine the Euler tour
               if obj.Ncc==1
                   [obj.Eul{1}, obj.De{1}, obj.H{1}, ~, obj.pa{1}] = GraphUtils.eulertour(obj.root, obj.ch{1});
               else
                   [obj.Eul{i}, obj.De{i}, obj.H{i}, ~, obj.pa{i}] = GraphUtils.eulertour(root_comp, obj.ch{i});
               end
               
               % Run the Range Minimum Query (RMQ) on the De matrix
               obj.I{i} = rmq(obj.De{i});
           end
       end  % end init
              
       % Send messages in queue; underlying graph: Gaussian MRF
       function send_g(obj,queue)
           for idx = 1:size(queue,1)
               src   = queue(idx,1);
               trg   = queue(idx,2);
               h_mrg = obj.h{src}      + sum(cat(3,obj.msg_h{src}{:}),3);
               J_mrg = obj.J{src,src}  + sum(cat(3,obj.msg_J{src}{:}),3);
               S     = -obj.J{trg,src}/(J_mrg - obj.msg_J{src}{obj.nghb2idx{src}(trg)});
               
               obj.msg_h{trg}{obj.nghb2idx{trg}(src)} = S*(h_mrg - obj.msg_h{src}{obj.nghb2idx{src}(trg)});
               obj.msg_J{trg}{obj.nghb2idx{trg}(src)} = S*obj.J{src,trg};
           end
       end  % end send_g
       
       % Send messages in queue; underlying graph: Discrete MRF (Sum-product)
       function send_d(obj,queue)
           for idx = 1:size(queue,1)
               src = queue(idx,1);
               trg = queue(idx,2);
               blf = obj.phi{src}.*prod(cell2mat(obj.msg{src}), 2);
               if iscell(obj.psi)
                   if obj.edgeOrig(src) ~= 0 && obj.edge2idx{obj.edgeOrig(src)}(trg) ~= 0
                       psi_cur = obj.psi{obj.edge2idx{obj.edgeOrig(src)}(trg)};
                   else
                       psi_cur = obj.psi{obj.edge2idx{obj.edgeOrig(trg)}(src)}';
                   end
               else
                   if obj.edgeOrig(src) ~= 0 && obj.edge2idx{obj.edgeOrig(src)}(trg) ~= 0
                       psi_cur = obj.psi;
                   else
                       psi_cur = obj.psi';
                   end
               end
               
               obj.msg{trg}{obj.nghb2idx{trg}(src)} = psi_cur'*(blf./obj.msg{src}{obj.nghb2idx{src}(trg)});  % here, parentheses DO matter! trust me!
               obj.msg{trg}{obj.nghb2idx{trg}(src)}(isnan(obj.msg{trg}{obj.nghb2idx{trg}(src)})) = 0;
               
               % Scale it up if it is too small
               if sum(obj.msg{trg}{obj.nghb2idx{trg}(src)}) < 1
                   obj.msg{trg}{obj.nghb2idx{trg}(src)} = obj.msg{trg}{obj.nghb2idx{trg}(src)}./sum(obj.msg{trg}{obj.nghb2idx{trg}(src)});
               end
           end
       end  % end send_d
       
       % Send messages in queue; underlying graph: Discrete MRF (Max-product)
       function send_m(obj,queue)
           for idx = 1:size(queue,1)
               src = queue(idx,1);
               trg = queue(idx,2);
               
               blf = log(obj.phi{src}) + sum(cell2mat(obj.msg{src}), 2);
               
               if iscell(obj.psi)
                   if obj.edgeOrig(src) ~= 0 && obj.edge2idx{obj.edgeOrig(src)}(trg) ~= 0
                       psi_cur = obj.psi{obj.edge2idx{obj.edgeOrig(src)}(trg)};
                   else
                       psi_cur = obj.psi{obj.edge2idx{obj.edgeOrig(trg)}(src)}';
                   end
               else
                   if obj.edgeOrig(src) ~= 0 && obj.edge2idx{obj.edgeOrig(src)}(trg) ~= 0
                       psi_cur = obj.psi;
                   else
                       psi_cur = obj.psi';
                   end
               end
               
               [max_val, max_idx] = max(bsxfun(@plus, log(psi_cur), blf - obj.msg{src}{obj.nghb2idx{src}(trg)}));  % subtract m_{trg->src}
               
               obj.msg{trg}{obj.nghb2idx{trg}(src)}     = max_val(:);  % m_{trg->src}
               obj.msg_idx{trg}{obj.nghb2idx{trg}(src)} = max_idx(:);  % indices of that maximize x_{src} for every value of x_{trg}
           end
       end  % end send_m
   end
   
   methods (Access = private, Static = true)
       % Determine the message schedule from src node to trg node
       function queue = schedule(lca, src, trg, pa)
           queue     = [];
           stack_src = [];
           stack_trg = [];
           
           node = src;
           while node ~= lca
               stack_src = [stack_src, node];
               node      = pa(node);
           end
           
           node = trg;
           while node ~= lca
               stack_trg = [stack_trg, node];
               node      = pa(node);
           end
           
           stack = [stack_src, lca, fliplr(stack_trg)];
           
           if length(stack) > 1
               queue = zeros(length(stack)-1,2);
               for i = 1:length(stack)-1
                   queue(i,:) = [stack(i)  stack(i+1)];
               end
           end
       end  % end schedule
   end
   
   methods
       % If nargin==3, varargin{1} = E, varargin{2} = phi, varargin{3} = psi
       % If nargin==2, varargin{1} = h, varargin{2} = J
       function obj = AdaBP(varargin)
           if any(strcmp(varargin,'max'))
               num_argin = nargin + 1;
           else
               num_argin = nargin;
           end
           
           if mod(num_argin,2) == 0
               obj.h = varargin{1};
               obj.J = varargin{2};
               obj.isGaussian = true;
               
               if ~iscell(obj.h),  obj.h = num2cell(obj.h);  end
               if ~iscell(obj.J),  obj.J = num2cell(obj.J);  end
               
               % Get number of nodes
               obj.N = length(obj.h);
               
               % Determine the edges
               tmp = cellfun(@(x) eq(x,0), obj.J, 'UniformOutput',false);
               tmp = cellfun(@all,         tmp  , 'UniformOutput',false);
               tmp = cellfun(@all,         tmp  , 'UniformOutput',false);
               tmp = ~cell2mat(tmp);
               [i,j] = find(triu(tmp~=0)-eye(size(tmp)));
               obj.E = [i,j];
               
           else
               obj.E   = varargin{1};
               obj.phi = varargin{2};
               obj.psi = varargin{3};
               obj.isGaussian = false;
               if iscell(obj.phi)   % if phi is cell, get information about node labels from phi
                   obj.N = length(obj.phi);
                   if any(strcmp(varargin,'lbl'))  % labels of the nodes in the graph
                       lbl = sort(varargin{find(strcmp(varargin,'lbl'))+1});
                   else
                       lbl = unique(obj.E(:));
                   end
                   
                   % Find the labels of nodes which participate in edges
                   lbl_in_E = unique(obj.E(:));
                   [isInGraph,new_lbl] = ismember(lbl_in_E, lbl);
                   
                   if ~all(isInGraph)
                       error('Not all nodes that participate in an edge appear in the graph. Check your labeling again.');
                   end
                   
                   % Relabel nodes from 1 to N
                   obj.E = GraphUtils.relabel(obj.E, new_lbl);
               else % if phi is common across all nodes, get information about node labels from E
                   obj.N   = length(unique(obj.E(:)));
                   tmp     = cell(1,obj.N);
                   tmp(:)  = {obj.phi};
                   obj.phi = tmp;                   
                   
                   % Relabel nodes from 1 to N
                   obj.E = GraphUtils.relabel(obj.E);
               end
           end
           
           % Determine the root
           if ~any(strcmp(varargin,'root'))
               obj.root = randi(obj.N);  % Choose a root randomly
           else
               obj.root = varargin{find(strcmp(varargin,'root'))+1};
           end
           
           % Determine whether marginals (Sum-Product) or MAP sequence is sought (Max-Product)
           obj.isSumProd = true;
           if ~obj.isGaussian  && any(strcmp(varargin,'max'))
               obj.isSumProd = false;
           end
           
           obj.init();
       end  % end AdaBP constructor
       
       % Delete all class variables
       function delete(obj)
           obj.N        = [];  obj.Ne = [];  obj.E    = [];  obj.root = [];
           obj.pa       = [];  obj.ch = [];  obj.nghb = [];
           obj.edge2idx = [];            obj.nghb2idx = [];
           obj.Eul      = [];  obj.De = [];     obj.H = [];  obj.I = [];
           obj.isGaussian = [];
           obj.phi = [];  obj.psi = [];  obj.nghb  = [];  obj.msg   = [];
           obj.h   = [];  obj.J   = [];  obj.msg_h = [];  obj.msg_J = [];
       end  % end delete
       
       % Update the node potential of specified node
       function update(obj, k, varargin)
           if obj.isGaussian
               Y_k   = varargin{1};
               C_k   = varargin{2};
               muW_k = varargin{3};
               R_k   = varargin{4};
               if isempty(muW_k),  muW_k = zeros(size(Y_k));  end
               tmp1       = C_k'/R_k;
               tmp2       = tmp1*C_k;
               tmp2       = (tmp2 + tmp2')/2;
               obj.h{k}   = obj.h{k} + tmp1*(Y_k - muW_k);
               obj.J{k,k} = obj.J{k,k} + tmp2;
               obj.J{k,k} = (obj.J{k,k} + obj.J{k,k}')/2;
           else
               obj.w_cur  = k;
               chi        = varargin{1};
               obj.phi{k} = obj.phi{k}.*chi;
               if sum(obj.phi{k}) < 1
                   obj.phi{k} = obj.phi{k}/sum(obj.phi{k});
               end
           end
           % Store this node as the last encountered measurement node from
           % this component
           obj.prev_w(obj.node2compID(k)) = obj.node2compNode(k);
       end  % end update
       
       % Set the node potential of specified node
       function setNodePot(obj, k, varargin)
           if obj.isGaussian
               if ~isempty(varargin{1})
                   obj.h{k} = varargin{1};
               end
               if length(varargin)==2 && ~isempty(varargin{2})
                   obj.J{k,k} = varargin{2};
                   obj.J{k,k} = (obj.J{k,k} + obj.J{k,k}')/2;
               end
           else
               obj.w_cur  = k;
               obj.phi{k} = varargin{1};
               if sum(obj.phi{k}) < 1
                   obj.phi{k} = obj.phi{k}/sum(obj.phi{k});
               end
           end
           % Store this node as the last encountered measurement node from
           % this component
           obj.prev_w(obj.node2compID(k)) = obj.node2compNode(k);
       end  % end setNodePot       
       
       % Propagate messages from src node to trg node
       % If we have a vector of source and target nodes, we do for every
       % pair of (source,target) nodes
       % is_ww: indicates whether the propagation concerns two measurement
       % nodes (ww), or a measurement and an inference node (wv)
       % If it is the former (is_ww=true), we still need to propagate messages 
       % from the last encountered measurement node of the component that 
       % the current measurement node belongs to.
       % If it is called with an output, it returns the number of messages 
       % that are propagated between src_vec and trg_vec.
       % If src_vec, trg_vec are scalars, the number of messages equals the
       % path length between these two nodes.
       function num_msg = propagate(obj, src_vec, trg_vec, is_ww)
           if nargin ~= 4, is_ww = false;  end
           
           if length(src_vec) == 1 && length(trg_vec) > 1  % we send from a meaurement node to multiple marginal nodes
               tmp = trg_vec;  trg_vec = src_vec;  src_vec = tmp;
               reverseFlag = true;
           else
               reverseFlag = false;
           end
           
           queue_all = cell(length(trg_vec),1);  % contains the messages that need to be sent
           cnt       = 1;                        % from the source nodes to each target node
           for trg = trg_vec
               for src = src_vec
                   compID_src = obj.node2compID(src);    % connected component of source node
                   compID_trg = obj.node2compID(trg);    % connected component of target node
                   
                   % Relative index of the source in its connected component
                     % if source node belongs to the same component as target node, this is the source node
                     if compID_src == compID_trg
                         compNode_src = obj.node2compNode(src);
                         % if source node does not belong to the same component as target node, 
                         % source node is the most recent encountered measurement node
                         % from the component of the target node
                     elseif is_ww && obj.prev_w(compID_trg) ~= 0
                         compNode_src = obj.prev_w(compID_trg);
                         % a) if we have a wv pair AND source node does not belong to
                         % the same component as target node, then we do not need to
                         % propagate.
                         % b) if we have a ww pair AND source node does not belong to
                         % the same component as target node AND prev_w(compID_trg)=0,
                         % this means it is the first time we enter the component of
                         % target node and hence no propagation is necessary.
                     else
                         obj.queue = [];
                         num_msg   = [];
                         return;
                     end
                     
                     % Relative index of the target in its connected component
                     compID       = compID_trg;
                     compNode_trg = obj.node2compNode(trg);
                     
                     if obj.H{compID}(compNode_src) < obj.H{compID}(compNode_trg)
                         low = obj.H{compID}(compNode_src);
                         upp = obj.H{compID}(compNode_trg);
                     else
                         low = obj.H{compID}(compNode_trg);
                         upp = obj.H{compID}(compNode_src);
                     end
                     
                     k = floor(log2(upp - low + 1));
                     if obj.De{compID}(obj.I{compID}(low,k+1)) <= obj.De{compID}(obj.I{compID}(upp-2^k+1,k+1))
                         lca_idx = obj.I{compID}(low,k+1);
                     else
                         lca_idx = obj.I{compID}(upp-2^k+1,k+1);
                     end
                     
                     % Determine the lca(el_src,el_trg)
                     lca = obj.Eul{compID}(lca_idx);
                     
                     % Determine the schedule from src to trg
                     queue_cur = AdaBP.schedule(lca, compNode_src, compNode_trg, obj.pa{compID});
                     
                     % If the graph is comprised of more than one connected
                     % components, node indices in the queue (for this component)
                     % should be switched back to their original indices in the graph
                     if obj.Ncc > 1
                         queue_cur = obj.sameCCNodes{compID}(queue_cur);
                     end
                     
                     queue_all{cnt} = [queue_all{cnt}; queue_cur];                     
               end
               cnt = cnt + 1;
           end
           
           % Omit repeating messages
           if length(src_vec) > 1 || length(trg_vec) > 1
               queue_tmp = [];
               for i = 1:length(trg_vec)
                   if reverseFlag
                       queue_tmp = [queue_tmp; fliplr(unique(flipud(queue_all{i}),'rows','stable'))];
                   else
                       queue_tmp = [queue_tmp; flipud(unique(flipud(queue_all{i}),'rows','stable'))];
                   end
               end
               obj.queue = queue_tmp;
           else
               obj.queue = queue_all{1};
           end
           
           if nargout > 0,  num_msg = max(0, size(obj.queue,1));  end
           
           % Send the messages from src to trg
           if obj.isSumProd
               if obj.isGaussian
                   obj.send_g(obj.queue);
               else
                   obj.send_d(obj.queue);
               end
           else
               obj.send_m(obj.queue);
           end
       end  % end propagate
       
       % Reset all messages
       function reset(obj)
           if obj.isGaussian
               for i = 1:obj.N
                   if ~isempty(obj.msg_h{i})
                       obj.msg_h{i}(:) = {zeros(size(obj.msg_h{i}{1}))};
                   end
                   
                   if ~isempty(obj.msg_J{i})
                       obj.msg_J{i}(:) = {zeros(size(obj.msg_J{i}{1}))};
                   end
               end
           else
               for i = 1:obj.N
                   if ~isempty(obj.msg{i})
                       obj.msg{i}(:) = {ones(size(obj.msg{i}{1}))};
                   end
               end
           end
       end  % end reset
       
       % Evaluate the marginal of specified node
       function varargout = eval_mrg(obj, k)
           if obj.isGaussian
               h_mrg     = obj.h{k}    + sum(cat(3,obj.msg_h{k}{:}),3);
               J_mrg     = obj.J{k,k}  + sum(cat(3,obj.msg_J{k}{:}),3);
               varargout = {h_mrg, J_mrg};
           else
               blf       = obj.phi{k}.*prod(cell2mat(obj.msg{k}), 2);
               blf       = blf./sum(blf);
               varargout = {blf};
           end
       end  % end eval_mrg
       
       % Evaluate the MAP sequence
       function map_seq = mapseq(obj)
           if obj.w_cur == 0
               map_seq = obj.seq_prv;
               return;
           end
           
           if length(obj.sameCCNodes) == 1
               b              = log(obj.phi{obj.w_cur}) + sum(cell2mat(obj.msg{obj.w_cur}),2);
               [~, mroot_idx] = max(b);
               
               if isempty(obj.queue)
                   isDirty = logical(sparse(obj.N,obj.N));
               else
                   isDirty = logical(sparse(obj.queue(:,1),obj.queue(:,2),1,obj.N,obj.N));
               end
               
               obj.seq_prv = btrack(obj.w_cur, mroot_idx, obj.nghb, obj.msg_idx, 'lbl2idx', obj.lbl2idx{1}, 'nghb2idx', obj.nghb2idx, 'seq_prv', obj.seq_prv, 'isDirty', isDirty);  % backtracking
               map_seq     = obj.seq_prv;
           else
               % Retrieve the component of the current walk element (all other components would remain unaffected)
               rootNode    = obj.w_cur;
               compID      = obj.node2compID(rootNode);
               phi_cc      = obj.phi(obj.sameCCNodes{compID});
               nghb_cc     = obj.nghb(obj.sameCCNodes{compID});
               msg_cc      = obj.msg(obj.sameCCNodes{compID});
               msg_idx_cc  = obj.msg_idx(obj.sameCCNodes{compID});
               nghb2idx_cc = obj.nghb2idx(obj.sameCCNodes{compID});
               seq_prv_cc  = obj.seq_prv(obj.sameCCNodes{compID});
               
               b = log(phi_cc{obj.node2compNode(rootNode)}) + sum(cell2mat(msg_cc{obj.node2compNode(rootNode)}),2);
               [~, mroot_idx] = max(b);
               
               if isempty(obj.queue)
                   isDirty = logical(sparse(obj.N,obj.N));
               else
                   isDirty = logical(sparse(obj.queue(:,1),obj.queue(:,2),1,obj.N,obj.N));
               end
               
               map_seq_cc = btrack(rootNode, mroot_idx, nghb_cc, msg_idx_cc, 'lbl2idx', obj.lbl2idx{compID}, 'nghb2idx', nghb2idx_cc, 'seq_prv', seq_prv_cc, 'isDirty', isDirty);  % backtracking
               obj.seq_prv(obj.sameCCNodes{compID}) = map_seq_cc;
               map_seq = obj.seq_prv;
           end
       end  % end mapseq
       
       function log_lik = get_loglik(obj,seq)
           log_lik = loglik(seq, obj.E, obj.phi, obj.psi);
       end  % end get_loglik
   end % methods
end % classdef