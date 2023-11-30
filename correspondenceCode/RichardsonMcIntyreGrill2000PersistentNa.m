function k = RichardsonMcIntyreGrill2000PersistentNa

% Kinetic scheme
% Na channel
% Richardson, McIntyre, Grill 2000 Med Biol Eng Comput
% 
% Resting potential -82mV, Temperature 20oC

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'Fast Na+';
k.chanfield =                               'nafast';
k.issodium =                                true;

k.cond.value =                              0.005;
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              50;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            1;
k.gates.temp =                              20;
k.gates.label =                             {'m'};
k.gates.numbereach =                        3;
k.gates.alpha.equ =                         {'0.186*(V+48.4)./(1-exp(-(V+48.4)./10.3))'};
k.gates.beta.equ =                          {'0.0086*(-(V+42.7))./(1-exp((V+42.7)./9.16))'};
k.gates.alpha.q10 =                         2.2;
k.gates.beta.q10 =                          2.2;