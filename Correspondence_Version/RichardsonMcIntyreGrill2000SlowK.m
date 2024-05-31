function k = RichardsonMcIntyreGrill2000SlowK

% Kinetic scheme
% K channel
% Richardson, McIntyre, Grill 2000 Med Biol Eng Comput
% 
% Resting potential -82mV, Temperature 20oC

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'K+';
k.chanfield =                               'k';
k.issodium =                                false;

k.cond.value =                              0.08;
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              -84;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            1;
k.gates.temp =                              20;
k.gates.label =                             {'s'};
k.gates.numbereach =                        1;
k.gates.alpha.equ =                         {'0.00122*(V+19.5)./(1-exp(-(V+19.5)/23.6))'};
k.gates.beta.equ =                          {'0.000739*(-(V+87.1))./(1-exp((V+87.1)/21.8))'};
k.gates.alpha.q10 =                         3;
k.gates.beta.q10 =                          3;