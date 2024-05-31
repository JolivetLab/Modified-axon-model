function k = Kv_dimer_JPN

% Kinetic scheme
% K channel
% M J Christie, J P Adelman, J Douglass, R A North
% 
% Resting potential -82mV, Temperature 24oC

k = createKineticStructure;

k.scheme =                                  'HH';
k.channame =                                'KvILT';
k.chanfield =                               'k';
k.issodium =                                false;

k.cond.value =                              0.0; %
k.cond.units =                              {2,' S','cm',[1 -2]};
k.erev.value =                              -65;
k.erev.units =                              {1, 'mV', 1};

k.vrest.value =                             -82;%-82;
k.vrest.units =                             {1, 'mV', 1};

k.gates.number =                            2;
k.gates.temp =                              22; 
k.gates.label =                             {'m','h'};
k.gates.numbereach =                        [4, 1];

k.gates.inf.equ =                           {'(1 + exp((-48-V)./6)).^(-1/2)', ...
                                             '(1-0.5).*((1+exp(-(-71-V)./10)).^(-1))+0.5'};
k.gates.tau.equ =                           {'2.9+(0.031.*exp((V+60)./6)+0.083.*exp(-(V + 60)./45)).^(-1)', ...
                                             '96.5+1000.*(0.52.*exp((V+60)./20)+0.52.*exp(-(V+60)./8)).^(-1)'};

k.gates.alpha.equ =                         {'((1+exp((-48-V)./6)).^(-1/2))./(2.9+(0.031.*exp((V+60)./6)+0.083.*exp(-(V+60)./45)).^(-1))',...
                                             '((1-0.5).*((1+exp(-(-71-V)./10)).^(-1))+0.5)./(96.5+1000.*(0.52.*exp((V+60)./20)+0.52.*exp(-(V+60)./8)).^(-1))'};
k.gates.beta.equ =                          {'1./(2.9+(0.031.*exp((V+60)./6)+0.083.*exp(-(V+60)./45)).^(-1))-((1+exp((-48-V)./6).^(-1/2))./(2.9+(0.031.*exp((V+60)./6)+0.083.*exp(-(V+60)./45)).^(-1)))',...
                                             '1./(96.5+1000.*(0.52.*exp((V+60)./20)+0.52.*exp(-(V+60)./8)).^(-1)) - ((1-0.5).*((1+exp(-(-71-V)./10)).^(-1))+0.5)./(96.5+1000.*(0.52.*exp((V+60)./20)+0.52.*exp(-(V+60)./8)).^(-1))'};

k.gates.alpha.q10 =                         [3, 3];
k.gates.beta.q10 =                          [3, 3];