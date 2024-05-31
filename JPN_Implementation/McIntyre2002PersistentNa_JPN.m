function k = McIntyre2002PersistentNa_JPN

% Kinetic scheme
% Na channel
% Richardson, McIntyre, Grill 2000 Med Biol Eng Comput
% 
% Resting potential -82mV, Temperature 20oC
% Further Notes: changed by Sonia Balan to fit parameters from McIntyre
% 2002 fast Na parameters

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'Fast Na+';
k.chanfield =                               'nafast';
k.issodium =                                true;

k.cond.value =                              0.005;
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              60;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -72;%-82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            1;
k.gates.temp =                              36;
k.gates.label =                             {'m'};
k.gates.numbereach =                        3;
k.gates.alpha.equ =                         {'0.0353*(V+27)./(1-exp(-(V+27)/10.2))'};
k.gates.beta.equ =                          {'0.000883*(-(V+34))./(1-exp((V+34)./10))'};
k.gates.alpha.q10 =                         2.2;
k.gates.beta.q10 =                          2.2;