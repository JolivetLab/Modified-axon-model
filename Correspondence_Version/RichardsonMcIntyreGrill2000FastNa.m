function k = RichardsonMcIntyreGrill2000FastNa

% Kinetic scheme
% Na channel
% Richardson, McIntyre, Grill 2000 Med Biol Eng Comput
% 
% Resting potential -82mV, Temperature 20oC

% Notes: very similar to Schwarz, Reid, Bostock 1995

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'Fast Na+';
k.chanfield =                               'nafast';
k.issodium =                                true;

k.cond.value =                              3;
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              50;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            2;
k.gates.temp =                              20;
k.gates.label =                             {'m','h'};
k.gates.numbereach =                        [3, 1];
k.gates.alpha.equ =                         {'1.86*(V+25.4)./(1-exp(-(V+25.4)./10.3))',...
                                             '0.0336*(-(V+118))./(1-exp((V+118)./11))'};
k.gates.beta.equ =                          {'0.086*(-(V+29.7))./(1-exp((V+29.7)./9.16))',...
                                             '2.3./(1+exp(-(V+35.8)./13.4))'};
k.gates.alpha.q10 =                         [2.2, 2.9];
k.gates.beta.q10 =                          [2.2, 2.9];