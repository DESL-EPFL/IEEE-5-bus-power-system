function [sl,pv,pq,Psp,Qsp,Vsp,Pmax,Qmin,Qmax,Vmin,Vmax,Imax,Y,YL,YT,kpv,kqv,kpf,kqf] = createGridFromFiles(busdata,linedata)
    %% Bus data
    % File format:
    %   1: Type: 1: slack, 2: PV, 3: PQ
    %   2,3: Psp,Qsp: Specified active, reactive powers [p.u.]
    %   4: Vsp: Specified voltage profile [p.u.]
    %   5: Pmax: Max gen power [p.u.]
    %   6,7: Qmin,Qmax: Min, max reactive powers [p.u.]
    %   8,9: Vmin,Vmax: Min, max voltages [p.u.]
    %   10,11: kpv,kqv: Pload = Psp(V/Vsp)^kpv, Qload = Qsp(V/Vsp)^kqv
    %   12,13: kpf,kqf: Pload = Psp(1+kpf*df), Qload = Qsp(1+kqf*df)  
    
    type = busdata(:,1);
    Psp = busdata(:,2);
    Qsp = busdata(:,3);
    Vsp = busdata(:,4);
    Pmax = busdata(:,5);
    Qmin = busdata(:,6);
    Qmax = busdata(:,7);
    Vmin = busdata(:,8);
    Vmax = busdata(:,9);
    kpv = busdata(:,10);
    kqv = busdata(:,11);
    kpf = busdata(:,12);
    kqf = busdata(:,13);
    
    sl = find(type == 1);
    pv = find(type == 2);
    pq = find(type == 3);
    
    n = length(sl) + length(pv) + length(pq);
    
    %% Line data
    % File format:
    %   1,2: fromBus, toBus
    %   3,4,5: r [p.u.], x [p.u.], b [p.u.]
    %   6: Imax [p.u.]
    
    fromBus = linedata(:,1);
    toBus = linedata(:,2);
    r = linedata(:,3);
    x = linedata(:,4)*50/60;
    b = linedata(:,5)*50/60;
    Imax_tmp = linedata(:,6);
    
    Imax = zeros(n,n);
    for i = 1:length(Imax_tmp)
        Imax(fromBus(i),toBus(i)) = Imax_tmp(i);
        Imax(toBus(i),fromBus(i)) = Imax_tmp(i);
    end
    
    [Y,YL,YT] = YMatrix(fromBus,toBus,r,x,b);
end

