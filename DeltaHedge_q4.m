function [pnlTimeBased, pnlMoveBased, transTimeBased, transMoveBased, Ntrades] = ...
            DeltaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot)

    % Asset price
    S = NaN(Nsims, Nsteps+1);
    S(:,1) = S0;
    transTimeBased = zeros(Nsims,1);
    transMoveBased = zeros(Nsims,1);

    
    %%%%% time-based
    BankTimeBased = NaN(Nsims, Nsteps+1);  % dollar value in the bank account 
    DeltaTimeBased = NaN(Nsims, Nsteps+1); % Delta position  
    NumSTimeBased = NaN(Nsims, Nsteps+1);  % position in the asset
    
    % initial positions
    DeltaTimeBased(:,1) = Delta(S0, T, pos, K,  rf, sigma, OptType);
    NumSTimeBased(:,1)  = -DeltaTimeBased(:,1);
    
    transTimeBased(:) = abs(NumSTimeBased(:,1)) .* trans;
    BankTimeBased(:,1)  = -Price(S0, T, pos, K, rf, sigma, OptType) ...
                            - NumSTimeBased(:,1).*S0 - transTimeBased;
                        
   
    %%%%% moved-based
    NumSMoveBased = NaN(Nsims, Nsteps+1); % position in asset
    KeepDelta = false(Nsims, Nsteps+1);   % stores whether to adjust Delta in move-based hedge 
    BankMoveBased = NaN(Nsims, Nsteps+1); % dollar value in the bank accounts
    DeltaBd = zeros(Nsims, 2); % Delta levels surrounding current delta hedge position

    % initial positions
    NumSMoveBased(:,1) = -DeltaTimeBased(:,1);
    
    transMoveBased(:) = abs(NumSTimeBased(:,1)) .* trans;
    BankMoveBased(:,1) = -Price(S0, T, pos, K,  rf, sigma, OptType) ...
                            - NumSMoveBased(:,1).*S0 - transMoveBased;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    % step size
    dt = T/Nsteps;
    t = [0:dt:T];
    
    for i = 2 : Nsteps
 
        % simulate new asset price (using P-measure since we want real scenarios)
        Z = randn(Nsims,1);
        S(:,i) = S(:,i-1).*exp( (mu-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z);

      
        %%%%%% time-based delta-hedging
        
        % compute new Deltas
        DeltaTimeBased(:,i) = Delta(S(:,i), T-t(i), pos, K,  rf, sigma, OptType);        
        NumSTimeBased(:,i) = -DeltaTimeBased(:,i);
           
        % update bank accounts
        BankTimeBased(:,i) =  BankTimeBased(:,i-1)*exp(rf*dt) ...
                            - ( NumSTimeBased(:,i) - NumSTimeBased(:,i-1) ) .* S(:,i) ...
                            - abs( NumSTimeBased(:,i-1) - NumSTimeBased(:,i) )  * trans;
        
        transTimeBased(:) = transTimeBased(:) + abs(NumSTimeBased(:,i-1) - NumSTimeBased(:,i)) * trans;   
                        
       
        %%%%%% moved-based delta-hedging
        
        % compute bounding delta-levels
        DeltaBd(:,1) = NumSMoveBased(:,i-1) - DDelta;
        DeltaBd(:,2) = DeltaBd(:,1) + 2 * DDelta; 
        
        % has the Delta moved beyond the bands 
        KeepDelta(:,i) =  ( NumSTimeBased(:,i) > DeltaBd(:,1)) ...
                          & ( NumSTimeBased(:,i) < DeltaBd(:,2));
        
        % update move-based positions if moved enough, otherwise keep old positions
        NumSMoveBased(:,i) =   NumSMoveBased(:, i-1) .* KeepDelta(:,i)...
                              + NumSTimeBased(:, i) .* (~KeepDelta(:,i));
                
        % update bank accounts
        BankMoveBased(:,i) =  BankMoveBased(:,i-1) * exp(rf*dt) ...
                            - (NumSMoveBased(:,i) - NumSMoveBased(:,i-1)) .* S(:,i) ...
                            - abs(NumSMoveBased(:,i-1) - NumSMoveBased(:,i)) * trans;
        
        transMoveBased(:) = transMoveBased(:) + abs(NumSMoveBased(:,i-1) - NumSMoveBased(:,i)) * trans;              
        
    end

    % treat last-period differently (no rebalancing, just unwind) 
    % and discount to current time
    
    % simulate new asset price
    Z = randn(Nsims,1);
    S(:,end) = S(:,end-1).*exp( (mu-0.5*sigma^2)*dt + sigma*sqrt(dt)*Z);


    pnlTimeBased ...
               =  ( BankTimeBased(:,Nsteps)*exp(rf*dt) ...
                   + NumSTimeBased(:,Nsteps).*S(:,end)...
                   - max(0,K - S(:,end)));
                                         
    pnlMoveBased ...
               =  ( BankMoveBased(:,Nsteps)*exp(rf*dt) ...
                   + NumSMoveBased(:,Nsteps).*S(:,end) ...
                   - max(0,K - S(:,end)));

  
   
    
    % plot results....
    
    if MakeSamplePathPlot == true
        
        % Time Based...
        fig = figure(100);
        clf(fig);
        
        subplot(2,1,1);
        plot(t, S(1:10,:));
        xlabel('time', 'Fontsize',24)
        ylabel('Spot Price', 'Fontsize',24)
        title('Time Based -- Spot price paths', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);
        
        subplot(2,1,2);
        plot(t, NumSTimeBased(1:10,:))
        xlabel('time', 'Fontsize',24)
        ylabel('Delta', 'Fontsize',24)
        title('Time Based -- Delta Position', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);
        
        % Move Based
        fig = figure(101);
        clf(fig);
        
        subplot(2,1,1);
        plot(t,S(1:10,:));
        xlabel('time', 'Fontsize',24)
        ylabel('Spot Price', 'Fontsize',24)
        title('Move Based -- Spot price paths', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);
        
        subplot(2,1,2);
        plot(t, NumSMoveBased(1:10,:))
        xlabel('time', 'Fontsize',24)
        ylabel('Delta', 'Fontsize',24)
        title('Move Based -- Delta Position', 'Fontsize',24)
        set(gca,'fontsize',18);
        xlim([0 T]);
    end
        
    % historgram of terminal PnL
    qtl1 = quantile(pnlTimeBased, 1);
    qtl2 = quantile(pnlTimeBased, 0.000);
    
    bins = [qtl2 : (qtl1-qtl2)/51: qtl1];
    
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
    
    
    qtl2 = quantile(pnlMoveBased-pnlTimeBased, 0.999);
    qtl1 = quantile(pnlMoveBased-pnlTimeBased, 0.001);
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
    fprintf('Delta Hedging! \n');
    fprintf('  P&L Time Based: mean = %8f    stdev = %8f...\n',[mean(pnlTimeBased) std(pnlTimeBased)]);
    fprintf('  P&L Move Based: mean = %8f    stdev = %8f...\n\n', [mean(pnlMoveBased) std(pnlMoveBased)]);
    fprintf('  Diff P&L : mean = %8f    stdev = %8f...\n\n', [mean(pnlMoveBased-pnlTimeBased) std(pnlMoveBased-pnlTimeBased)/sqrt(Nsims)]);
    
    
    
    % compare replication with payoff
    terminalReplTimeBased = BankTimeBased(:,Nsteps)*exp(rf*dt) + NumSTimeBased(:,Nsteps).*S(:,end);
    terminalReplMoveBased = BankMoveBased(:,Nsteps)*exp(rf*dt) + NumSMoveBased(:,Nsteps).*S(:,end);
    
    
    Sup = quantile(S(:,end), 0.999);
    Slow = quantile(S(:,end), 0.001);
    Srep = [Slow : (Sup-Slow)/101 : Sup]';
    
    fig = figure(4);
    clf(fig);
    
    scatter(S(:,Nsteps+1), terminalReplTimeBased, '.r')
    hold
    plot( Srep, -Price(Srep, 1e-5, pos, K, rf, sigma, OptType), '-k','LineWidth',1);
    hold
    title('Time Based Replication','fontsize',24);
    xlabel('$S_T$','interpreter','latex','fontsize',18);
    ylabel('$D_{t_{n-1}}e^{r\Delta t_n}+\alpha_{t_{n-1}}S_{t_n}$','interpreter','latex','fontsize',18);
    set(gca,'fontsize',18);
    
    fig = figure(5);
    clf(fig);
    
    scatter(S(:,Nsteps+1), terminalReplMoveBased, '.r')
    hold
    plot( Srep, -Price(Srep, 1e-5, pos, K, rf, sigma, OptType), '-k','LineWidth',1);
    hold
    title('Move Based Replication','fontsize',24);
    xlabel('$S_T$','interpreter','latex','fontsize',18);
    ylabel('$D_{t_{n-1}}e^{r\Delta t_n}+\alpha_{t_{n-1}}S_{t_n}$','interpreter','latex','fontsize',18);
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


% price, delta
function result = Delta(S, T, pos, K, rf, sigma, OptType)

    result = zeros(length(S),1);

    for n = 1 : length(K)
        dp = (log(S./K) + (rf + 1/2*sigma.^2).*T) ./ (sigma.*sqrt(T));
        delta = normcdf(dp, 0, 1) - 1;
        result = result + pos(n)*delta;
     
    end

end

function result = Price(S, T, pos, K, rf, sigma, OptType)

    result = zeros(length(S),1);

    for n = 1 : length(K)
        dp = ( log(S./K) + (rf + 1/2*sigma.^2).*T) ./ (sigma.*sqrt(T));
        dm = dp - sigma.*sqrt(T);
        price =  - S .* normcdf(-dp, 0, 1) + exp(-rf*T) .* K .* normcdf(-dm, 0, 1);
        result = result + pos(n)*price;

    end

end

