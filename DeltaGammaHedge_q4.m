function [pnlTimeBased, pnlMoveBased, transTimeBased, transMoveBased, Ntrades] = DeltaGammaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, trans_e, trans_op, DDelta, Nsims, Nsteps, MakeSamplePathPlot)

    % option maturity for hedging option: using at-the-money
    Thdg = T + 0.25;
    Shdg = S0;

    % Asset price
    S = NaN(Nsims, Nsteps+1);
    S(:,1) = S0;
    
    transTimeBased = zeros(Nsims,1);
    transMoveBased = zeros(Nsims,1);
   

    
    %%%%% time-based
    BankTimeBased = NaN(Nsims, Nsteps+1);  % dollar value in the bank account 
    DeltaTimeBased = NaN(Nsims, Nsteps+1); % Delta position  
    GammaTimeBased = NaN(Nsims, Nsteps+1); % Delta position  
    NumSTimeBased = NaN(Nsims, Nsteps+1);  % position in the asset
    NumOptTimeBased = NaN(Nsims, Nsteps+1);  % position in hedging option
    
    % initial positions
    
    DeltaTimeBased(:,1) = -abs(Delta(S0, T, pos, K,  rf, sigma, OptType)); %delta of put, negative
    Delta_call = abs(Delta(S0, Thdg, 1, Shdg,  rf, sigma, true)); %delta of call, positive
    GammaTimeBased(:,1) = abs(Gamma(S0, T, pos, K,  rf, sigma, OptType)); %gamma of put, positive
    Gamma_call = abs(Gamma(S0, Thdg, 1, Shdg, rf, sigma, true)); %gamma of call, positive
    
    NumOptTimeBased(:,1) = GammaTimeBased(:,1)./ Gamma_call;    
    NumSTimeBased(:,1) = DeltaTimeBased(:,1) - NumOptTimeBased(:,1) .* Delta_call;

    Price_put = abs(Price(S0, T, pos, K,   rf, sigma, OptType));
    Price_call = abs(Price(S0, Thdg, 1, Shdg,   rf, sigma, true));
    
    BankTimeBased(:,1) = Price_put ...
                        - NumSTimeBased(:,1).* S0 ...
                        - trans_e * abs(NumSTimeBased(:,1))...
                        - NumOptTimeBased(:,1) .* Price_call ...
                        - trans_op * abs(NumOptTimeBased(:,1));
    
    transTimeBased(:) = trans_e * abs(NumSTimeBased(:,1)) + trans_op * abs(NumOptTimeBased(:,1));

    
    %%%%% moved-based
    NumSMoveBased = NaN(Nsims, Nsteps+1); % position in asset
    NumOptMoveBased = NaN(Nsims, Nsteps+1); % position in hedging option
    KeepDelta = false(Nsims, Nsteps+1);   % stores whether to adjust Delta in move-based hedge 
    BankMoveBased = NaN(Nsims, Nsteps+1); % dollar value in the bank accounts
    DeltaBd = zeros(Nsims, 2); % Delta levels surrounding current delta hedge position

    % initial positions
    NumSMoveBased(:,1) = NumSTimeBased(:,1);
    NumOptMoveBased(:,1) = NumOptTimeBased(:,1);
    BankMoveBased(:,1) = Price_put ...
                        - NumSMoveBased(:,1).* S0 ...
                        - abs(NumSMoveBased(:,1))*trans_e ...
                        - NumOptMoveBased(:,1) .* Price_call ...
                        - abs(NumOptMoveBased(:,1)) * trans_op;
                    
    transMoveBased(:) = abs(NumSMoveBased(:,1))*trans_e + abs(NumOptMoveBased(:,1)) * trans_op;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % step size
    dt = T/Nsteps;
    t = [0:dt:T];
    
    for i = 2 : Nsteps
 
        % simulate new asset price (using P-measure since we want real scenarios)
        Z = randn(Nsims,1);
        S(:,i) = S(:,i-1).*exp( (mu-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z);

       
        %%%%%% time-based delta-hedging
        
        % compute new Gammas
        GammaTimeBased(:,i) = abs(Gamma(S(:,i), T-t(i), pos, K, rf, sigma, OptType));
        Gamma_call = abs(Gamma(S(:,i), Thdg-t(i), 1, Shdg,  rf, sigma, true));
        NumOptTimeBased(:,i) =  GammaTimeBased(:,i) ./ Gamma_call;
        
        % compute new Deltas
        DeltaTimeBased(:,i) = -abs(Delta(S(:,i), T-t(i), pos, K,  rf, sigma, OptType));
        Delta_call = abs(Delta(S(:,i), Thdg-t(i), 1, Shdg,  rf, sigma, true));
        NumSTimeBased(:,i) = DeltaTimeBased(:,i) - (NumOptTimeBased(:,i)) .* Delta_call;
                                  
        % update bank accounts
        Price_call = abs(Price(S(:,i), Thdg-t(i), 1, Shdg,  rf, sigma, true));
        BankTimeBased(:,i) =  BankTimeBased(:,i-1)*exp(rf*dt) ...
                            - ( NumSTimeBased(:,i) - NumSTimeBased(:,i-1) ) .* S(:,i) ...
                            - abs( NumSTimeBased(:,i-1) - NumSTimeBased(:,i) ) * trans_e ...
                            - ( NumOptTimeBased(:,i) - NumOptTimeBased(:,i-1) ) .* Price_call ...
                            - abs( NumOptTimeBased(:,i-1) - NumOptTimeBased(:,i) ) * trans_op; 
        
        transTimeBased(:) = transTimeBased(:) ...
                            + abs( NumSTimeBased(:,i-1) - NumSTimeBased(:,i) ) * trans_e ...
                            + abs( NumOptTimeBased(:,i-1) - NumOptTimeBased(:,i) ) * trans_op;
               

        %%%%%% moved-based delta-hedging
        
        % compute bounding delta-levels
        DeltaBd(:,1) = NumSMoveBased(:,i-1) - DDelta;
        DeltaBd(:,2) = DeltaBd(:,1) + 2 * DDelta; 
        
        % has the Delta moved beyond the bands ?
        KeepDelta(:,i) =    ( NumSTimeBased(:,i) > DeltaBd(:,1) ) ...
                          & ( NumSTimeBased(:,i) < DeltaBd(:,2));
        
        % update move-based positions if moved enough, otherwise keep old positions
        NumSMoveBased(:,i) =   NumSMoveBased(:, i-1) .* KeepDelta(:,i)...
                              + NumSTimeBased(:, i) .* (~KeepDelta(:,i));
                

        NumOptMoveBased(:,i) =   NumOptMoveBased(:, i-1) .* KeepDelta(:,i)...
                              + NumOptTimeBased(:, i) .* (~KeepDelta(:,i));
                          
                          
        % update bank accounts
        BankMoveBased(:,i) =  BankMoveBased(:,i-1)*exp(rf*dt) ...
                            - ( NumSMoveBased(:,i) - NumSMoveBased(:,i-1) ) .* S(:,i) ...
                            - abs( NumSMoveBased(:,i-1) - NumSMoveBased(:,i) ) * trans_e ...
                            - ( NumOptMoveBased(:,i) - NumOptMoveBased(:,i-1) ) .* Price_call ...                            
                            - abs( NumOptMoveBased(:,i-1) - NumOptMoveBased(:,i) ) * trans_op;
                        
        transMoveBased(:) = transMoveBased(:) ...
                            + abs( NumSMoveBased(:,i-1) - NumSMoveBased(:,i) ) * trans_e ...
                            + abs( NumOptMoveBased(:,i-1) - NumOptMoveBased(:,i) ) * trans_op;
        
    end

    % treat last-period differently (no rebalancing, just unwind) 
    % and discount to current time
    
    % simulate new asset price
    Z = randn(Nsims,1);
    S(:,end) = S(:,end-1).*exp( (mu-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z);
    
    BankTimeBased(1,end-1)*exp(rf*dt);
   
    terminalReplTimeBased =   BankTimeBased(:,end-1)*exp(rf*dt) ...
                            + NumSTimeBased(:,end-1).*S(:,end) ...
                            + NumOptTimeBased(:,end-1).*Price(S(:,end), (Thdg-T)+1e-5, 1, Shdg, rf, sigma, true);

    terminalReplMoveBased = BankMoveBased(:,end-1)*exp(rf*dt) ...
                            + NumSMoveBased(:,end-1).*S(:,end) ...
                            + NumOptMoveBased(:,end-1).*Price(S(:,end), (Thdg-T)+1e-5, 1, Shdg, rf, sigma, true);
    
    
    
    pnlTimeBased =  terminalReplTimeBased - max(0,K - S(:,end));
    pnlMoveBased =  terminalReplMoveBased - max(0,K - S(:,end));
    % plot results....
    
    if MakeSamplePathPlot
        
        % Time Based...
        fig = figure(100);
        clf(fig);
        
        subplot(3,1,1);
        plot(t,S(1:10,:));
        xlabel('time', 'Fontsize',24)
        ylabel('Spot Price', 'Fontsize',24)
        title('Time Based -- Spot price paths', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);
        
        subplot(3,1,2);
        plot(t,NumSTimeBased(1:10,:));
        xlabel('time', 'Fontsize',24)
        ylabel('Delta', 'Fontsize',24)
        title('Time Based -- Delta Position', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);
        
        subplot(3,1,3);
        plot(t,NumOptTimeBased(1:10,:))
        xlabel('time', 'Fontsize',24)
        ylabel('Gamma', 'Fontsize',24)
        title('Time Based -- Gamma Position', 'Fontsize',24)
        set(gca,'fontsize',18);        
        xlim([0 T]);
        
        % Move Based
        fig = figure(101);
        clf(fig);
        
        subplot(3,1,1);
        plot(t,S(1:10,:));
        xlabel('time', 'Fontsize',24)
        ylabel('Spot Price', 'Fontsize',24)
        title('Move Based -- Spot price paths', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);

        subplot(3,1,2);
        plot(t,NumSMoveBased(1:10,:))
        xlabel('time', 'Fontsize',24)
        ylabel('Delta', 'Fontsize',24)
        title('Move Based -- Delta Position', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);

        subplot(3,1,3);
        plot(t,NumOptMoveBased(1:10,:));
        xlabel('time', 'Fontsize',24)
        ylabel('Gamma', 'Fontsize',24)
        title('Time Based -- Gamma Position', 'Fontsize',24)
        set(gca,'fontsize',18);        
        xlim([0 T]);
        
        
    end
        
    % historgram of terminal PnL
    qtl1 = quantile(pnlTimeBased, 0.995);
    qtl2 = quantile(pnlTimeBased, 0.005);
    
    qtl = max(abs([qtl1 qtl2]));
    bins = [-qtl : 2*qtl/51: qtl];
    
    fig = figure(1);
    clf(fig);
    
    [ freqTime binsTime ] = hist(pnlTimeBased, bins);
    hist(pnlTimeBased, bins);
    title('Time Based P&L','fontsize',24);
    xlim([bins(1) bins(end)]);
    set(gca,'fontsize',18);
    
    fig = figure(201);
    clf(fig); 
    
    [ freqMove binsMove ] = hist(pnlMoveBased, bins);
    hist(pnlMoveBased, bins);
    title('Move Based P&L','fontsize',24);
    xlim([bins(1) bins(end)]);
    set(gca,'fontsize',18);

    
    fig = figure(3);
    clf(fig);
    
    harray=bar(binsTime, [freqTime' freqMove'],'group');
    set(harray(2),'faceColor','r','barwidth',2);
    legend('Time Base','Move Based');
    title('Compare P&Ls','fontsize',24);
    xlim([bins(1) bins(end)]);
    set(gca,'fontsize',18);
    
    
    
    qtl2 = quantile(pnlMoveBased-pnlTimeBased, 0.995);
    qtl1 = quantile(pnlMoveBased-pnlTimeBased, 0.005);
    bins = [qtl1 : (qtl2-qtl1)/51: qtl2];

    fig = figure(10);
    %[ freqMove binsMove ] = 
    hist(pnlMoveBased-pnlTimeBased, bins);
    %hist(pnlMoveBased, bins);
    title('$\Delta$P\&L','fontsize',24,'interpreter','latex');
    xlim([bins(1) bins(end)]);
    set(gca,'fontsize',18);

    
    
    
    % compare P and L
    fprintf('\n');
    fprintf('Delta Gamma Hedging!\n');
    fprintf('  P&L Time Based: mean = %8f    stdev = %8f...\n',[mean(pnlTimeBased) std(pnlTimeBased)]);
    fprintf('  P&L Move Based: mean = %8f    stdev = %8f...\n\n', [mean(pnlMoveBased) std(pnlMoveBased)]);
    fprintf('  Diff P&L : mean = %8f    stdev = %8f...\n\n', [mean(pnlMoveBased-pnlTimeBased) std(pnlMoveBased-pnlTimeBased)]);
   
    
    % compare replication with payoff    
    Sup = quantile(S(:,end), 0.995);
    Slow = quantile(S(:,end), 0.005);
    Srep = [Slow : (Sup-Slow)/101 : Sup]';
    
    fig = figure(4);
    clf(fig);
    
    scatter(S(:,Nsteps+1), terminalReplTimeBased, '.r')
    hold
    plot( Srep, -Price(Srep, 1e-5, pos, K, rf, sigma, OptType), '-k','LineWidth',1);
    hold
    title('Time Based Replication','fontsize',24);
    xlabel('$S_T$','interpreter','latex','fontsize',18);
    ylabel('$D_{t_{n-1}}e^{r\Delta t_n}+\alpha_{t_{n-1}}S_{t_n}+\beta_{t_{n-1}}g_{t_n}$','interpreter','latex','fontsize',18);
    set(gca,'fontsize',18);
    
    fig = figure(5);
    clf(fig);
    
    scatter(S(:,Nsteps+1), terminalReplMoveBased, '.r')
    hold
    plot( Srep, -Price(Srep, 1e-5, pos, K, rf, sigma, OptType), '-k','LineWidth',1);
    hold
    title('Move Based Replication','fontsize',24);
    xlabel('$S_T$','interpreter','latex','fontsize',18);
    ylabel('$D_{t_{n-1}}e^{r\Delta t_n}+\alpha_{t_{n-1}}S_{t_n}+\beta_{t_{n-1}}g_{t_n}$','interpreter','latex','fontsize',18);
    set(gca,'fontsize',18);

    
    % plot histogram of number of trades using moved-based trading
    fig = figure(6);
    clf(fig);

    Ntrades = sum(1-KeepDelta,2);
    hist(Ntrades, [0:1:Nsteps]); xlim([0 Nsteps]);
    title('Histogram of # Trades','fontsize',24);
    set(gca,'fontsize',18);
    
    
    
    fprintf(' Mean Time Based Trans Cost = %8f   stdev = %8f \n',[mean(transTimeBased) std(transTimeBased)]);
    fprintf(' Mean Move Based Trans Cost = %8f   stdev = %8f \n', [mean(transMoveBased) std(transMoveBased)]);
    fprintf(' Mean Number of Transactions = %8f \n', mean(Ntrades));

    
end


% price, delta, gamma
function result = Gamma(S, T, pos, K, rf, sigma, OptType)

    result = zeros(length(S),1);

    for n = 1 : length(K)
        if OptType(n)
            result = result + pos(n)*CallGamma(S, K(n), T, rf, sigma);
        else
            result = result + pos(n)*PutGamma(S, K(n), T, rf, sigma);
        end
    end

end


function result = Delta(S, T, pos, K, rf, sigma, OptType)

    result = zeros(length(S),1);

    for n = 1 : length(K)
        if OptType(n)
            result = result + pos(n)*CallDelta(S, K(n), T, rf, sigma);
        else
            result = result + pos(n)*PutDelta(S, K(n), T, rf, sigma);
        end
    end

end

function result = Price(S, T, pos, K, rf, sigma, OptType)

    result = zeros(length(S),1);

    for n = 1 : length(K)
        if OptType(n)
            result = result + pos(n)*CallPrice(S, K(n), T, rf, sigma);
        else
            result = result + pos(n)*PutPrice(S, K(n), T, rf, sigma);
        end
    end

end

function price = CallPrice(S, K, T, rf, sigma)

    dp = (log(S./K) + (rf + 1/2*sigma.^2)*T) ./ (sigma*sqrt(T));
    dm = dp - sigma*sqrt(T);
    
    price = S .* normcdf(dp, 0, 1) - K .* exp(-rf*T) .* normcdf(dm, 0, 1);
    
end

function delta = CallDelta(S, K, T, rf, sigma)

    dp = ( log(S ./ K) + (rf + 1/2 * sigma.^2).* T ) ./ (sigma.*sqrt(T));
    delta = normcdf(dp, 0, 1);
    
end


function price = PutPrice(S, K, T, rf, sigma)

    dp = ( log(S./K) + (rf + 1/2*sigma.^2).*T) ./ (sigma.*sqrt(T));
    dm = dp - sigma.*sqrt(T);
    
    price =  - S .* normcdf(-dp, 0, 1) + exp(-rf*T) .* K .* normcdf(-dm, 0, 1);
    
end

function delta = PutDelta(S, K, T, rf, sigma)

    dp = (log(S./K) + (rf + 1/2*sigma.^2).*T) ./ (sigma.*sqrt(T));
    delta = normcdf(dp, 0, 1) - 1;
    
end

function gamma = PutGamma(S, K, T, rf, sigma)

    gamma = CallGamma(S, K, T, rf, sigma);
    
end

function gamma = CallGamma(S, K, T, rf, sigma)

    dp = ( log(S./K) + (rf + 1/2*sigma.^2).*T) ./ (sigma.*sqrt(T));
    gamma = normpdf(dp, 0, 1) ./ (sigma .* S .* sqrt(T) );
    
end