%% Parameter Setup

S0 = 100;
mu = 0.1;
rf = 0.02;
sigma = 0.2;
T = 1/4;

trans = 0.005;
trans_op = 0.01;
DDelta = 0.05/2;
Nsims = 10000;
Nsteps = 100;
MakeSamplePathPlot = false;


%The put option to be hedged
pos= [-1];
K =[100];
OptType=[false];

%% Delta Hedge

MakeSamplePathPlot = false;

[pnlTimeBased_d, pnlMoveBased_d] = ...
    DeltaHedge(S0, mu, rf, sigma, T, pos, K, OptType, ...
               trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot);


%% Delta-Gamma Hedge

MakeSamplePathPlot = false;

[pnlTimeBased_dg, pnlMoveBased_dg] = ...
    DeltaGammaHedge(S0, mu, rf, sigma, T, pos, K, OptType, ...
                    trans, trans_op, DDelta, Nsims, Nsteps, MakeSamplePathPlot);                                        

%% Real volatility vs. risk-neutral volatility - Delta hedge

MakeSamplePathPlot = false;
sigma_p = [0.1 0.15 0.2 0.25 0.3];


[pnlTimeBased_q3d, pnlMoveBased_q3d] = ...
    DeltaHedge_q3(S0, mu, rf, sigma, sigma_p, T, pos, K, OptType, ...
                  trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot);

%% Real volatility vs. risk-neutral volatility - Delta-Gamma hedge

MakeSamplePathPlot = false;
sigma_p = [0.1 0.15 0.2 0.25 0.3];

[pnlTimeBased_q3dg, pnlMoveBased_q3dg] = ...
    DeltaGammaHedge_q3(S0, mu, rf, sigma, sigma_p, T, pos, K, OptType, ...
                       trans, trans_op, DDelta, Nsims, Nsteps, MakeSamplePathPlot);

%% Change in Rebalancing band

DDelta = 0.01/2;

[pnlTimeBased_d4a, pnlMoveBased_d4a, transTimeBased_d4a, transMoveBased_d4a, Ntrades_d4a] ...
    = DeltaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, ...
               trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot);
          
%%           
DDelta = 0.05/2;
[pnlTimeBased_d4b, pnlMoveBased_d4b, transTimeBased_d4b, transMoveBased_d4b, Ntrades_d4b] ...
    = DeltaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, ...
               trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot);

%%
DDelta = 0.1/2;
[pnlTimeBased_d4c, pnlMoveBased_d4c, transTimeBased_d4c, transMoveBased_d4c, Ntrades_d4c] ...
    = DeltaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, ...
               trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot);

%%         
DDelta = 0.5/2;
[pnlTimeBased_d4d, pnlMoveBased_d4d, transTimeBased_d4d, transMoveBased_d4d, Ntrades_d4d] ...
    = DeltaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, ...
               trans, DDelta, Nsims, Nsteps, MakeSamplePathPlot);

%% Delta-Gamma
DDelta = 0.01/2;

[pnlTimeBased_dg4a, pnlMoveBased_dg4a, transTimeBased_dg4a, transMoveBased_dg4a, Ntrades_dg4a] ...
    = DeltaGammaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, trans, trans_op, DDelta, Nsims, Nsteps, MakeSamplePathPlot);
                                        
%%
DDelta = 0.05/2;
[pnlTimeBased_dg4b, pnlMoveBased_dg4b, transTimeBased_dg4b, transMoveBased_dg4b, Ntrades_dg4b] ...
    = DeltaGammaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, trans, trans_op, DDelta, Nsims, Nsteps, MakeSamplePathPlot);
          
%%
DDelta = 0.1/2;
[pnlTimeBased_dg4c, pnlMoveBased_dg4c, transTimeBased_dg4c, transMoveBased_dg4c, Ntrades_dg4c] ...
    = DeltaGammaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, trans, trans_op, DDelta, Nsims, Nsteps, MakeSamplePathPlot);
          
%%
DDelta = 0.5/2;
[pnlTimeBased_dg4d, pnlMoveBased_dg4d, transTimeBased_dg4d, transMoveBased_dg4d, Ntrades_dg4d] ...
    = DeltaGammaHedge_q4(S0, mu, rf, sigma, T, pos, K, OptType, trans, trans_op, DDelta, Nsims, Nsteps, MakeSamplePathPlot);
 
