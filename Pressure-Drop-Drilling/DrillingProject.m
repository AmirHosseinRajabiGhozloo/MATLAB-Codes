clc;
clear;
close all;
   teta300 = input(['enter your teta 300 = ' ]);
   teta600 = input(['enter your teta 600 = ' ]);
   rho = input(['enter your density = ' ]);
    qm = input(['enter your maximum deby = ']);
    dip = input(['enter your in diameter of drill pipe = ' ]);
    dop = input(['enter your out diameter of drill pipe = ' ]);
    dic = input(['enter your in diameter of drill collar = ' ]);
    doc = input(['enter your out diameter of drill collar = ' ]);
    dc = input(['enter your casing in diameter = ' ]);
    dh = input(['enter your open hole diameter = ' ]);
    db = input(['enter your bit diameter = ' ]);
    ldp = input(['enter your drill pipe length = ']); %please pay attantion to that ldp+ldc=lc+hl
    ldc = input(['enter your collar pipe length = ']);
    lc = input(['enter your casing length = ']);
    hl = input(['enter your open hole lenght = ']);
    dn1 = input(['enter your first nozzle diameter = ']);
    dn2 = input(['enter your second nozzle diameter = ']);
    dn3 = input(['enter your third nozzle diameter = ']);
    cd = input (['enter your cd = ']);
    ft = input('what is your fluid type (newtonian / bingham)= '); %please give your answer between two quotation
    switch ft
        case 'newtonian'
for q = 1:1:qm
    vip = PipeVelocity(q,dip);
    NReip = NewtonianReynoldsPipe(rho,vip,dip,teta300);
    if NReip < 2100
        ippl = LaminarNewtonianPressureLossPipe(vip,dip,teta300);
    else
        ippl = NewtonianTurbulentPressureLossPipe(rho,vip,dip,teta300); 
    end
    vic = PipeVelocity(q,dic);
    NReic = NewtonianReynoldsPipe(rho,vic,dic,teta300);
    if NReic < 2100
        icpl = LaminarNewtonianPressureLossPipe(vic,dic,teta300);
    else
        icpl = NewtonianTurbulentPressureLossPipe(rho,vic,dic,teta300);
    end
    bpl = BitPressureLoss(rho,q,cd,dn1,dn2,dn3);
   if  ldc >= hl
           vac1 = AnnulusVelocity(q,dh,doc);
           NReac1 = NewtonianReynoldsAnnulus(rho,vac1,doc,dh,teta300);
           if NReac1 < 2100
               ocpl1 = LaminarNewtonianPressureLossAnnulus(vac1,dh,doc,teta300);
           else
               ocpl1 = NewtonianTurbulentPressureLossAnnulus(rho,vac1,doc,dh,teta300);
           end
           vac2 = AnnulusVelocity(q,dc,doc);
           NReac2 = NewtonianReynoldsAnnulus(rho,vac2,doc,dc,teta300);
           if NReac2 < 2100
               ocpl2 = LaminarNewtonianPressureLossAnnulus(vac2,dc,doc,teta300);
           else
               ocpl2 = NewtonianTurbulentPressureLossAnnulus(rho,vac2,doc,dc,teta300);
           end
           vap = AnnulusVelocity(q,dc,dop);
           NReap = NewtonianReynoldsAnnulus(rho,vap,dop,dc,teta300);
           if NReap < 2100
               oppl = LaminarNewtonianPressureLossAnnulus(vap,dc,dop,teta300);
           else
               oppl = NewtonianTurbulentPressureLossAnnulus(rho,vap,dop,dc,teta300);
           end
           tpl(1,q) = bpl+(ocpl1*hl)+(ocpl2*(ldc-hl))+(oppl*ldp)+(ippl*ldp)+(icpl*ldc);
           fpl(1,q) = (ocpl1*hl)+(ocpl2*(ldc-hl))+(oppl*ldp)+(ippl*ldp)+(icpl*ldc);
           Bpl(1,q) = bpl ;
       else
           vac = AnnulusVelocity(q,dh,doc);
           NReac = NewtonianReynoldsAnnulus(rho,vac,doc,dh,teta300);
           if NReac < 2100
               ocpl = LaminarNewtonianPressureLossAnnulus(vac,dh,doc,teta300);
           else
               ocpl = NewtonianTurbulentPressureLossAnnulus(rho,vac,doc,dh,teta300);
           end
           vap1 = AnnulusVelocity(q,dh,dop);
           NReap1 = NewtonianReynoldsAnnulus(rho,vap1,dop,dh,teta300);
           if NReap1 < 2100
               oppl1 = LaminarNewtonianPressureLossAnnulus(vap1,dh,dop,teta300);
           else
               oppl1 = NewtonianTurbulentPressureLossAnnulus(rho,vap1,dop,dh,teta300);
           end
           vap2 = AnnulusVelocity(q,dc,dop);
           NReap2 = NewtonianReynoldsAnnulus(rho,vap2,dop,dc,teta300);
           if NReap2 < 2100
               oppl2 = LaminarNewtonianPressureLossAnnulus(vap2,dc,dop,teta300);
           else
               oppl2 = NewtonianTurbulentPressureLossAnnulus(rho,vap2,dop,dc,teta300);
           end
           tpl(1,q) = bpl+(ocpl*ldc)+(oppl1*(hl-ldc))+(oppl2*lc)+(ippl*ldp)+(icpl*ldc);
           fpl(1,q) = (ocpl*ldc)+(oppl1*(hl-ldc))+(oppl2*lc)+(ippl*ldp)+(icpl*ldc);
           Bpl(1,q) = bpl;
   end
end
q=1:1:qm;
        figure(1);
        plot(q,tpl);
        xlabel('q')
        ylabel('delta P')
        figure(2);
        plot(q,fpl);
        xlabel('q')
        ylabel('delta P')
        figure(3);
        plot(q,Bpl);
        xlabel('q')
        ylabel('delta P')
        
        figure(4);
        plot(q,tpl);
        xlabel('q')
        ylabel('delta P')
        hold on;
        plot(q,fpl);
        plot(q,Bpl); 
        case 'bingham'
for q = 1:1:qm
    vip = PipeVelocity(q,dip);
    NReip = binghamReynoldsPipe(rho,vip,dip,teta300,teta600);
    if NReip < 2100
        ippl = LaminarBinghamPressureLossPipe(dip,vip,teta300,teta600);
    else
        ippl = BinghamTurbulentPressureLossPipe(rho,vip,teta300,teta600,dip);
    end
    vic = PipeVelocity(q,dic);
    NReic = binghamReynoldsPipe(rho,vic,dic,teta600,teta300);
    if NReic < 2100
        icpl = LaminarBinghamPressureLossPipe(dic,vic,teta300,teta600);
    else
        icpl = BinghamTurbulentPressureLossPipe(rho,vic,teta300,teta600,dic);
    end
    bpl = BitPressureLoss(rho,q,cd,dn1,dn2,dn3);
   if  ldc >= hl
           vac1 = AnnulusVelocity(q,dh,doc);
           NReac1 = BinghamReynoldsAnnulus(rho,vac1,dh,doc,teta600,teta300);
           if NReac1 < 2100
               ocpl1 = LaminarBinghamPressureLossAnnulus(doc,dh,vac1,teta300,teta600);
           else
               ocpl1 = BinghamTurbulentPressureLossAnnulus(rho,vac1,teta300,teta600,dh,doc);
           end
           vac2 = AnnulusVelocity(q,dc,doc);
           NReac2 = BinghamReynoldsAnnulus(rho,vac2,dc,doc,teta600,teta300);
           if NReac2 < 2100
               ocpl2 = LaminarBinghamPressureLossAnnulus(doc,dc,vac2,teta300,teta600);
           else
               ocpl2 = BinghamTurbulentPressureLossAnnulus(rho,vac2,teta300,teta600,dc,doc);
           end
           vap = AnnulusVelocity(q,dc,dop);
           NReap = BinghamReynoldsAnnulus(rho,vap,dc,dop,teta600,teta300);
           if NReap < 2100
               oppl = LaminarBinghamPressureLossAnnulus(dop,dc,vap,teta300,teta600);
           else
               oppl = BinghamTurbulentPressureLossAnnulus(rho,vap,teta300,teta600,dc,dop);
           end
           tpl(1,q) = bpl+(ocpl1*hl)+(ocpl2*(ldc-hl))+(oppl*ldp)+(ippl*ldp)+(icpl*ldc);
           fpl(1,q) = (ocpl1*hl)+(ocpl2*(ldc-hl))+(oppl*ldp)+(ippl*ldp)+(icpl*ldc);
           Bpl(1,q) = bpl ;
       else
           vac = AnnulusVelocity(q,dh,doc);
           NReac = BinghamReynoldsAnnulus(rho,vac,dh,doc,teta600,teta300);
           if NReac < 2100
               ocpl = LaminarBinghamPressureLossAnnulus(doc,dh,vac,teta300,teta600);
           else
               ocpl = BinghamTurbulentPressureLossAnnulus(rho,vac,teta300,teta600,dh,doc);
           end
           vap1 = AnnulusVelocity(q,dh,dop);
           NReap1 = BinghamReynoldsAnnulus(rho,vap1,dh,dop,teta600,teta300);
           if NReap1 < 2100
               oppl1 = LaminarBinghamPressureLossAnnulus(dop,dh,vap1,teta300,teta600);
           else
               oppl1 = BinghamTurbulentPressureLossAnnulus(rho,vap1,teta300,teta600,dh,dop);
           end
           vap2 = AnnulusVelocity(q,dc,dop);
           NReap2 = BinghamReynoldsAnnulus(rho,vap2,dc,dop,teta600,teta300);
           if NReap2 < 2100
               oppl2 = LaminarBinghamPressureLossAnnulus(dop,dc,vap2,teta300,teta600);
           else
               oppl2 = BinghamTurbulentPressureLossAnnulus(rho,vap2,teta300,teta600,dc,dop);
           end
           tpl(1,q) = bpl+(ocpl*ldc)+(oppl1*(hl-ldc))+(oppl2*lc)+(ippl*ldp)+(icpl*ldc);
           fpl(1,q) = (ocpl*ldc)+(oppl1*(hl-ldc))+(oppl2*lc)+(ippl*ldp)+(icpl*ldc);
           Bpl(1,q) = bpl;
   end
end
q=1:1:qm;
        figure(1);
        plot(q,tpl);
        xlabel('q')
        ylabel('delta P')
        figure(2);
        plot(q,fpl);
        xlabel('q')
        ylabel('delta P')
        figure(3);
        plot(q,Bpl);
        xlabel('q')
        ylabel('delta P')
        
        figure(4);
        plot(q,tpl);
        xlabel('q')
        ylabel('delta P')
        hold on;
        plot(q,fpl);
        plot(q,Bpl);
    end
      