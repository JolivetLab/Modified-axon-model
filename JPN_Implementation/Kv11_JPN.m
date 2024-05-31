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

k.cond.value =                              2.0; %
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              -65;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -72;%-82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            2;
k.gates.temp =                              24; 
k.gates.label =                             {'l','r'};
k.gates.numbereach =                        [1, 2];

k.gates.inf.equ =                           {'1.0./(1 + exp(-(V + 30.5)./11.3943))', ...
                                             '1.0./(1 + exp((V + 30.0)./27.3943))'};
k.gates.tau.equ =                           {'30.0./(1 + exp((V + 76.56)./26.1479))', ...
                                             '15000.0./(1 + exp((V + 160.56)/-100.0))'};

k.gates.alpha.equ =                         {'(1.0./(1 + exp(-(V + 30.5)./11.3943)))./(30.0./(1 + exp((V + 76.56)./26.1479)))',...
                                             '(1.0./(1 + exp((V + 30.0)./27.3943)))./(15000.0./(1 + exp((V + 160.56)/-100.0)))'};
k.gates.beta.equ =                          {'1./(30.0./(1 + exp((V + 76.56)./26.1479))) - (1.0./(1 + exp(-(V + 30.5)./11.3943)))./(30.0./(1 + exp((V + 76.56)./26.1479)))',...
                                             '1./(15000.0./(1 + exp((V + 160.56)/-100.0))) - (1.0./(1 + exp((V + 30.0)./27.3943)))./(15000.0./(1 + exp((V + 160.56)/-100.0)))'};

k.gates.alpha.q10 =                         [3, 3];
k.gates.beta.q10 =                          [3, 3];