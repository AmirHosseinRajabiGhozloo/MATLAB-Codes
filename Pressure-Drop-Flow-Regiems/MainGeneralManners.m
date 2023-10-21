clc;
clear;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ CHOOSING ONE OF GENERAL MANNERS ++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ 1 - HOMOGENEOUS +++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ 2 - LOCKHART _ MARTINELLI ++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ 3 - HAGEDORN & BROWN +++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

aaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
bbb = ('===>>  It is a General calculation method (not related to the well condition).  <<===');
ccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(aaa);
disp(bbb);
disp(ccc);

k = ['which kind of general manners do you prefer to calculate the pressure drop in two phase flow ?'];
disp(k);
n = input('[[[ 1 - Homogeneous   or   2 - Lockhart-Martinelli   or   3 - Hagedorn & Brown ]]] ===>  ');

switch n
    
    case 1
        Homogeneous_General
    case 2
        LockhartMartinelli_General
    case 3 
        HagedornBrown_General
    
end