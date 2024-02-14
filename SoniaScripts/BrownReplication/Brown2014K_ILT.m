function k = Brown2014K

% Kinetic scheme
% K channel
% Richardson, McIntyre, Grill 2000 Med Biol Eng Comput
% 
% Resting potential -82mV, Temperature 20oC
% Further Notes: changed by Sonia Balan to fit parameters from McIntyre
% 2002 fast Na parameters

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'K+';
k.chanfield =                               'k';
k.issodium =                                false;

k.cond.value =                              0.002;
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              -80;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -72;%-82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            1;
k.gates.temp =                              36;
k.gates.label =                             {'s'};
k.gates.numbereach =                        1;
k.gates.alpha.equ =                         {'0.3./(1+exp((V+53)./-5))'};
k.gates.beta.equ =                          {'0.03./(1+exp((V+90)./-1))'};
k.gates.alpha.q10 =                         3;
k.gates.beta.q10 =                          3;