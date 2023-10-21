clc;
clear;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ ALL WELLS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ DISTINGUISHING FLOW PATTERN AND PRESSURE DROP ++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

aaaa = ('Please enter the angle of the well from horizontal in degree. if you do not know please enter ===> 200');
disp(aaaa);
nn = input('Theta ====>>>  ');

if nn >= 0  &&  nn <= 5
    
    aaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    bbb = ('===>>  It is a Horizontal well or pipe.  <<===');
    ccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp(aaa);
    disp(bbb);
    disp(ccc);
    
elseif nn > 5  &&  nn < 80
    
    FlowPatternAndPressureDrop_Deviated;
    
elseif nn > 80  && nn <= 90
    
    FlowPatternAndPressureDrop_Vertical;
    
elseif nn == 200
    
    MainGeneralManners;
    
end
    