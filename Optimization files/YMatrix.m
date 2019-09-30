function [Y,YL,YT] = YMatrix(fromBus,toBus,r,x,b)
    n_nodes = max(max(fromBus),max(toBus));
    n_lines = length(r);
    
    YL = zeros(n_nodes);
    YT = zeros(n_nodes);
    
    for i = 1:n_lines
        YL(fromBus(i),toBus(i)) = 1/(r(i)+1j*x(i));
        YL(toBus(i),fromBus(i)) = 1/(r(i)+1j*x(i)); 
        YT(fromBus(i),toBus(i)) = 1j*b(i)/2;
        YT(toBus(i),fromBus(i)) = 1j*b(i)/2; 
    end
    
    Y = -YL + diag(YL*ones(n_nodes,1)) + diag(YT*ones(n_nodes,1));
end

