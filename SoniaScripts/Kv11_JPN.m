function k = Kv11_JPN

% Kinetic scheme
% K channel
% M J Christie, J P Adelman, J Douglass, R A North
% 
% Resting potential -82mV, Temperature 24oC

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'Kv1.1';
k.chanfield =                               'k';
k.issodium =                                false;

k.cond.value =                              0.002; % from Brown & Hamann
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              -65;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -72;%-82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            2;
k.gates.temp =                              36;
k.gates.label =                             {'l','r'};
k.gates.numbereach =                        [1, 2];
k.gates.alpha.equ =                         {'5.0./(1 + exp(-(V + 30.5)./11.3943))', ...
                                             '5.0./(1 + exp((V + 30.0)./27.3943))'};
k.gates.beta.equ =                          {'150.0./(1 + exp((V + 76.56)./26.1479))', ...
                                             '75000.0./(1 + exp((V + 160.56)/-100.0))'};
k.gates.alpha.q10 =                         [3, 3];
k.gates.beta.q10 =                          [3, 3];