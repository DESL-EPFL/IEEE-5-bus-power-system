function [Cp_V,Cq_V,Cp_I,Cq_I] = sensCoeffs(Y,YL,YT,V,I,Psp,Qsp,Vsp,kpv,kqv,sl,pv,pq,wrt_nodes)
% Compute sensitivity coefficients for |V| and |I| with respect to
% nodes specified in 'wrt_nodes'

    % Remove slacks
    V_no_sl = V;
    V_no_sl(sl) = [];
    Y_new = Y;
    Y_new(sl,:) = [];
    Y_new_no_sl = Y_new;
    Y_new_no_sl(:,sl) = [];
    G_new = real(Y_new);
    B_new = imag(Y_new);
    G_new_no_sl = real(Y_new_no_sl);
    B_new_no_sl = imag(Y_new_no_sl);

    Psp_no_sl = Psp;
    Qsp_no_sl = Qsp;
    Psp_no_sl(sl) = [];
    Qsp_no_sl(sl) = [];
    Vsp_no_sl = Vsp;
    kpv_no_sl = kpv;
    kqv_no_sl = kqv;
    Vsp_no_sl(sl) = [];
    kpv_no_sl(sl) = [];
    kqv_no_sl(sl) = [];
    
    nsl = length(sl);
    npv = length(pv);
    npq = length(pq);
    n = nsl+npv+npq;
    
    Cp_V = zeros(n-nsl,length(wrt_nodes)); % d|V|/dP
    Cq_V = zeros(n-nsl,length(wrt_nodes)); % d|V|/dQ
    Cp_V_com = zeros(n-nsl,length(wrt_nodes)); 
    Cq_V_com = zeros(n-nsl,length(wrt_nodes));
    
    %% Build system
    A_RR = (real(V_no_sl).*G_new_no_sl + imag(V_no_sl).*B_new_no_sl); % Real part of eq. <-> dV_R/dP
    A_RI = (-real(V_no_sl).*B_new_no_sl + imag(V_no_sl).*G_new_no_sl); % Real part of eq. <-> dV_I/dP
    A_IR = (imag(V_no_sl).*G_new_no_sl - real(V_no_sl).*B_new_no_sl); % Imaginary part of eq. <-> dV_R/dP
    A_II = (-imag(V_no_sl).*B_new_no_sl - real(V_no_sl).*G_new_no_sl); % Imaginary part of eq. <-> dV_I/dP

    % Diagonal terms
    A_RR = A_RR + diag(G_new*real(V)-B_new*imag(V)-Psp_no_sl./Vsp_no_sl.^kpv_no_sl.*kpv_no_sl.*abs(V_no_sl).^(kpv_no_sl-2).*real(V_no_sl));
    A_RI = A_RI + diag(G_new*imag(V)+B_new*real(V)-Psp_no_sl./Vsp_no_sl.^kpv_no_sl.*kpv_no_sl.*abs(V_no_sl).^(kpv_no_sl-2).*imag(V_no_sl));
    A_IR = A_IR + diag(-G_new*imag(V)-B_new*real(V)-Qsp_no_sl./Vsp_no_sl.^kqv_no_sl.*kqv_no_sl.*abs(V_no_sl).^(kqv_no_sl-2).*real(V_no_sl));
    A_II = A_II + diag(G_new*real(V)-B_new*imag(V)-Qsp_no_sl./Vsp_no_sl.^kqv_no_sl.*kqv_no_sl.*abs(V_no_sl).^(kqv_no_sl-2).*imag(V_no_sl));
    
    A = [A_RR A_RI;
         A_IR A_II];
     
    for k = 1:length(wrt_nodes)
        %% RHS
        j = wrt_nodes(k);
        
        bp_R = zeros(n-nsl,1); % Real part of eq. for P
        bp_I = zeros(n-nsl,1); % Imaginary part of eq. for P
        bp_R(k) = (abs(V(j))/Vsp(j))^kpv(j);
        
        bq_R = zeros(n-nsl,1); % Real part of eq. for Q
        bq_I = zeros(n-nsl,1); % Imaginary part of eq. for Q
        bq_I(k) = (abs(V(j))/Vsp(j))^kqv(j);
        
        %% Solve for P
        b = [bp_R;
             bp_I];
        
        Cp_V_tmp = A\b;
        Cp_V_com(:,k) = Cp_V_tmp(1:n-nsl) + 1j*Cp_V_tmp(n-nsl+1:end);
        
        Cp_V(:,k) = 1./abs(V_no_sl).*(real(V_no_sl).*real(Cp_V_com(:,k))+imag(V_no_sl).*imag(Cp_V_com(:,k)));
        
        %% Solve for Q
        b = [bq_R;
             bq_I];
        
        Cq_V_tmp = A\b;
        Cq_V_com(:,k) = Cq_V_tmp(1:n-nsl) + 1j*Cq_V_tmp(n-nsl+1:end);
        
        Cq_V(:,k) = 1./abs(V_no_sl).*(real(V_no_sl).*real(Cq_V_com(:,k))+imag(V_no_sl).*imag(Cq_V_com(:,k)));
    end
    
    %% Currents sensitivity coefficients
    Cp_I = zeros(n,n,length(wrt_nodes)); % d|I|/dP
    Cq_I = zeros(n,n,length(wrt_nodes)); % d|I|/dQ
    Cp_I_com = zeros(n,n,length(wrt_nodes));
    Cq_I_com = zeros(n,n,length(wrt_nodes));
    
    % Add slacks
    Cp_V_com_sl = zeros(n,length(wrt_nodes));
    Cq_V_com_sl = zeros(n,length(wrt_nodes));
    Cp_V_com_sl(sort([pv;pq]),:) = Cp_V_com;
    Cq_V_com_sl(sort([pv;pq]),:) = Cq_V_com;
    
    for k = 1:length(wrt_nodes)
        for i = 1:n
            for j = 1:n
                Cp_I_com(i,j,k) = YL(i,j)*(Cp_V_com_sl(i,k)-Cp_V_com_sl(j,k))+YT(i,j)*Cp_V_com_sl(i,k);
                Cq_I_com(i,j,k) = YL(i,j)*(Cq_V_com_sl(i,k)-Cq_V_com_sl(j,k))+YT(i,j)*Cq_V_com_sl(i,k);
                
                if i ~= j && abs(Y(i,j)) ~= 0
                    Cp_I(i,j,k) = 1/abs(I(i,j))*(real(I(i,j))*real(Cp_I_com(i,j,k))+imag(I(i,j))*imag(Cp_I_com(i,j,k)));
                    Cq_I(i,j,k) = 1/abs(I(i,j))*(real(I(i,j))*real(Cq_I_com(i,j,k))+imag(I(i,j))*imag(Cq_I_com(i,j,k)));
                end
            end
        end   
    end
end