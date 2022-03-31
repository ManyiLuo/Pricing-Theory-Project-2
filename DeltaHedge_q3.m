function [pnlTimeBased, pnlMoveBased] = DeltaHedge_q3(S0, mu, rf, sigma, sigma_p, T, pos, K, OptType, trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot)
    
    pnlTimeBased = [];
    pnlMoveBased = [];
    for j = 1:size(sigma_p,2)
        
        % Asset price
        S = NaN(Nsims, Nsteps+1);
        S(:,1) = S0;

     
        %%%%% time-based
        BankTimeBased = NaN(Nsims, Nsteps+1);  % dollar value in the bank account 
        DeltaTimeBased = NaN(Nsims, Nsteps+1); % Delta position  
        NumSTimeBased = NaN(Nsims, Nsteps+1);  % position in the asset

        % initial positions
        DeltaTimeBased(:,1) = Delta(S0, T, pos, K,  rf, sigma, OptType);
        NumSTimeBased(:,1)  = -DeltaTimeBased(:,1);
        BankTimeBased(:,1)  = -Price(S0, T, pos, K, rf, sigma, OptType) ...
                                - NumSTimeBased(:,1).*S0 - abs(NumSTimeBased(:,1)) * trans;
       
                            
        %%%%% moved-based
        NumSMoveBased = NaN(Nsims, Nsteps+1); % position in asset
        KeepDelta = false(Nsims, Nsteps+1);   % stores whether to adjust Delta in move-based hedge 
        BankMoveBased = NaN(Nsims, Nsteps+1); % dollar value in the bank accounts
        DeltaBd = zeros(Nsims, 2); % Delta levels surrounding current delta hedge position

        % initial positions
        NumSMoveBased(:,1) = -DeltaTimeBased(:,1);
        BankMoveBased(:,1) = -Price(S0, T, pos, K,  rf, sigma, OptType) ...
                                - NumSMoveBased(:,1).*S0 - abs(NumSTimeBased(:,1)) * trans;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % step size
        dt = T/Nsteps;
        t = [0:dt:T];

        for i = 2 : Nsteps

            % simulate new asset price (using P-measure since we want real scenarios)
            Z = randn(Nsims,1);
            S(:,i) = S(:,i-1).*exp( (mu-0.5*sigma_p(j)^2)*dt + sigma_p(j)*sqrt(dt)*Z);

       
            %%%%%% time-based delta-hedging

            % compute new Deltas
            DeltaTimeBased(:,i) = Delta(S(:,i), T-t(i), pos, K,  rf, sigma, OptType);        
            NumSTimeBased(:,i) = -DeltaTimeBased(:,i);

            % update bank accounts
            BankTimeBased(:,i) =  BankTimeBased(:,i-1)*exp(rf*dt) ...
                                - ( NumSTimeBased(:,i) - NumSTimeBased(:,i-1) ) .* S(:,i) ...
                                - abs( NumSTimeBased(:,i-1) - NumSTimeBased(:,i) )  * trans;

          
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


        end

        % treat last-period differently (no rebalancing, just unwind) 
        % and discount to current time

        % simulate new asset price
        Z = randn(Nsims,1);
        S(:,end) = S(:,end-1).*exp( (mu-0.5*sigma_p(j)^2)*dt + sigma_p(j)*sqrt(dt)*Z);


      
        pnlTimeBased(:,end+1) ...
                   =  ( BankTimeBased(:,Nsteps)*exp(rf*dt) ...
                       + NumSTimeBased(:,Nsteps).*S(:,end)...
                       - max(0,K - S(:,end)));

        pnlMoveBased(:,end+1) ...
                   =  ( BankMoveBased(:,Nsteps)*exp(rf*dt) ...
                       + NumSMoveBased(:,Nsteps).*S(:,end) ...
                       - max(0,K - S(:,end)));

        % compare P and L
        fprintf('Realized Volatility = %6.4f \n', sigma_p(j));
        fprintf('  P&L Time Based: mean = %8f    stdev = %8f...\n',[mean(pnlTimeBased(:,end)) std(pnlTimeBased(:,end))]);
        fprintf('  P&L Move Based: mean = %8f    stdev = %8f...\n\n', [mean(pnlMoveBased(:,end)) std(pnlMoveBased(:,end))]);
        fprintf('  Diff P&L : mean = %8f    stdev = %8f...\n\n', [mean(pnlMoveBased(:,end)-pnlTimeBased(:,end)) std(pnlMoveBased(:,end)-pnlTimeBased(:,end))]);
   
    end


    
    % KDE of PnL time based
    fig = figure(1);
    clf(fig);
    for j = 1: size(pnlTimeBased,2)

       [ freqTime x ] = ksdensity(pnlTimeBased(:,j));
       plot(x, freqTime,'LineWidth',3);
       hold on;
       title('KDE of Time Based P&L, Delta Hedging','fontsize',24);
       set(gca,'fontsize',18);
       xlim([-8 4]);
    end
    
    hold off;
    legend("sigma=10%","sigma=15%","sigma=20%","sigma=25%","sigma=30%",'Location','northwest');
    xlabel('P&L', 'Fontsize',18)
    ylabel('KDE', 'Fontsize',18)
    
    % KDE of PnL move based
    fig = figure(2);
    clf(fig);
    for j = 1: size(pnlMoveBased,2)

       [ freqTime x ] = ksdensity(pnlMoveBased(:,j));
       plot(x, freqTime,'LineWidth',3);
       hold on;
       title('KDE of Move Based P&L, Delta Hedging','fontsize',24);
       set(gca,'fontsize',18);
       xlim([-8 4]);
    end
    
    hold off;
    legend("sigma=10%","sigma=15%","sigma=20%","sigma=25%","sigma=30%",'Location','northwest');
    xlabel('P&L', 'Fontsize',18)
    ylabel('KDE', 'Fontsize',18)
    

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

