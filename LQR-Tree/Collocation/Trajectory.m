classdef Trajectory
    properties
        x_d
        u_d
        N
        dt
        S
        rho
    end
    
    methods
        function self = Trajectory(varargin)
            self.S = varargin{1};
            if ~isa(self.S, 'function_handle')
                self.S = @(t) self.S;
            end
            self.rho = varargin{2};
            self.N = varargin{3};
            self.dt = varargin{4};
            self.x_d = varargin{5};
            if ~isa(self.x_d, 'function_handle')
                self.x_d = @(t) self.x_d;
            end
            self.u_d = varargin{6};
          
            
        end
  
        
        
        function obj = TreeNode(modelHandle, attemptPosition, Q, R, TreeHead,u_max,N, dt) 
            %Everything after R should be excluded for root node
            %   if this node is not the root of tree, log parent node
            obj.position = attemptedPosition;
            obj.TimeStepSize = dt;
            obj.modelList = [];
            obj.modelList(1) = modelHandle(attemptedPosition, Q,R);
            if exist('TreeHead','var')
                % if not head, then trajectory to nearest node
                [closestNode, ~] = FindClosestNodeEuclideanRecursive(obj.position,TreeHead);

                [obj.TrajectoryX,obj.TrajectoryU] = collocate_trajectory(obj.modelList(1).dynamics(), obj.position, closestNode.position, u_max, N,dt);

                % instantiate models at each point along x trajectory
            else
                %if this node is head, do single position lqr
                obj.TimeStepSize = 0;
                
            end
        
            
        end

    end
    
    methods (Static)
        function [TreeNode, distance] = FindClosestNodeEuclideanRecursive(position, currentTreeNode)
            closest = currentTreeNode;
            closestDist = norm(currentTreeNode.position - position);
            for model = children
                [node, dist] = FindClosestNodeRecursive(position, model);
                if(dist<closestDist)
                    closest = node;
                    closestDist = dist;
                end
            end
            TreeNode = closest;
            distance = closestDist;
        end
        
        
    end
end


