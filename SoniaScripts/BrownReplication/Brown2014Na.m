function k = Brown2014Na

% Kinetic scheme
% Na channel
% Richardson, McIntyre, Grill 2000 Med Biol Eng Comput
% 
% Resting potential -82mV, Temperature 20oC

% Notes: very similar to Schwarz, Reid, Bostock 1995
% Further Notes: changed by Sonia Balan to fit parameters from McIntyre
% 2002 fast Na parameters

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'Fast Na+';
k.chanfield =                               'nafast';
k.issodium =                                true;

k.cond.value =                              0.05;
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              50;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -72;%-82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            2;
k.gates.temp =                              36;
k.gates.label =                             {'m','h'};
k.gates.numbereach =                        [3, 1];
k.gates.alpha.equ =                         {'6.57*(V+20.4)./(1-exp(-(V+20.4)./10.3))', ...
                                             '0.34*(-(V+114))./(1-exp((V+114)./11))'};
k.gates.beta.equ =                          {'0.304*(-(V+25.7))./(1-exp((V+25.7)./9.16))', ...
                                             '12.6./(1+exp(-(V+31.8)./13.4))'};
k.gates.alpha.q10 =                         [2.2, 2.9];
k.gates.beta.q10 =                          [2.2, 2.9];