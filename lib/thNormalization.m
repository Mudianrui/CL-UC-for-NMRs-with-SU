%[-pi,pi)
function thPi=thNormalization(th)
thPi = th-floor((th+pi)/(pi*2))*(pi*2);
