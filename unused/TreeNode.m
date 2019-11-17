classdef TreeNode

    properties
        position
        parent %single TreeNode instance
        children %vector of children TreeNode instances
        TimeStepSize %time step for trajectory
        TrajectoryX %trajectories of x
        TrajectoryU %trajectories of u
        modelList %list same size at trajectories with corresponding models with tvlqr
    end
    
    methods
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
                % I don't think collocation is setup to return x
                % correctly.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                [obj.TrajectoryX,obj.TrajectoryU] = collocate_trajectory(obj.modelList(1).dynamics(), obj.position, closestNode.position, u_max, N,dt);
                %!!!!!!!need to instantiate models around each x
                %point!!!!!!!!!!
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


