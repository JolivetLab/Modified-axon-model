function k = Kv12_JPN

% Kinetic scheme
% K channel
% M J Christie, J P Adelman, J Douglass, R A North
% 
% Resting potential -82mV, Temperature 24oC

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'Kv1.2';
k.chanfield =                               'k';
k.issodium =                                false;

k.cond.value =                              0.02; %
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              -65;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -72;%-82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            2;
k.gates.temp =                              20; 
k.gates.label =                             {'m','h'};
k.gates.numbereach =                        [1, 1];

k.gates.inf.equ =                           {'1.0./(1+ exp((V + 21.0)./-11.3943))', ...
                                             '1.0./(1+ exp((V + 22.0)./11.3943))'};
k.gates.tau.equ =                           {'150.0./(1+ exp((V + 67.56)./34.1479)) ', ...
                                             '15000.0./(1+ exp((V + 46.56)./-44.1479)) '};

k.gates.alpha.equ =                         {'(1.0./(1+ exp((V + 21.0)./-11.3943)))./(150.0./(1+ exp((V + 67.56)./34.1479)))',...
                                             '(1.0./(1+ exp((V + 22.0)./11.3943)))./(15000.0./(1+ exp((V +46.56)./-44.1479)))'};
k.gates.beta.equ =                          {'1./(150.0./(1+ exp((V + 67.56)./34.1479))) - (1.0./(1+ exp((V + 21.0)./-11.3943)))./(150.0./(1+ exp((V + 67.56)./34.1479)))',...
                                             '1./(15000.0./(1+ exp((V + 46.56)./-44.1479))) - (1.0./(1+ exp((V + 22.0)./11.3943)))./(15000.0./(1+ exp((V +46.56)./-44.1479)))'};

k.gates.alpha.q10 =                         [3, 3];
k.gates.beta.q10 =                          [3, 3];