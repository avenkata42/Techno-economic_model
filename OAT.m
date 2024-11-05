% function [ff]=SVSA_mod(bc,lb,ub)
% SVSAandprofit(bc,lb,ub)
% [ff]=SVSAandprofit(bc,lb,ub)
%%%%% above lines to integrate this code MS excel. 

% This code is a representative example to demonstrate how the "Single variable sensitivity analysis" (SVSA) or the "One-at-a-time" (OAT) is carried
% out in this study. The actual code is integrated with Aspen Plus for inputs on the component flow rates. Here sample calculations are used
% to approximate these numbers.

% clear all
% close all

% Code to calculate the Minimum selling price of Ethylene with and without
% selling of the other products

% Definition of sensitivity parameters
    titles=['Capture efficiency (%)', 'Capital cost (%)', "Threshold CO3 concentration (M)", "Membrane permeance (L/m^2/h/bar)[GO]", "CO3 Rejection (%)[GO]", 'Membrane cost ($/m^2)[GO]', "Membrane Lifetime (years)[GO]", 'Stack cost (%)', 'Single pass conversion (%)', 'Current density (mA/cm^2)', 'Price of Product (%)', 'Electricity cost ($/kWh)', 'Faradaic efficiency C2H4 (%)', 'Faradaic efficiency C2H5OH (%)',  'Cell Voltage (%)', 'Electrode area per cell (m^2)', 'H2 permeance (GPU)[HF]', 'C2H4 permeance (GPU)[HF]', 'Membrane cost ($/m^2)[HF]', 'Membrane lifetime (years)[HF]', 'Water flux (LMH)[Pervap]', 'Membrane selectivity [Pervap]', 'Membrane cost ($/m^2)[Pervap]', 'Membrane lifetime (years)[Pervap]'];
    
% Defining lower and upper bounds 
    lb=[-0.0625	-0.25	-0.25	-0.25	-0.333333333	-0.175	-0.666666667	-0.25	-0.428571429	-0.4	-0.15	-0.25	-0.32	-0.35	-0.15	-0.571428571	0	-0.333333333	-0.2	-0.083333333	-0.333333333	-0.15	-0.142857143	-0.25];
    ub=[0.125	0.15	0	0.25	0.055555556	0	0.666666667	0.1	0.857142857	1	0.15	3	0.4	0.25	0.43	0.714285714	0.363636364	0	0.4	0.166666667	0	0	0.071428571	0.25];

% Defining the base case values for the sensitivity parameters for the 3 cases 
    %     bc =[0.9	0.75	1.5	3.5	0.95	165	7	0.75	0.65	1	1.15	0.015	0.65 0.25 0.85	6	1500	0.8	20	7	1.5	10000	60 7]; %optimistic
    bc=[0.8	1	2	2.8	0.9	200	3	1	0.35	0.5	1	0.02	0.5	0.2	0.975	3.5	1100	1.2	25	6 1.5	10000	70	4];  %base case
    % bc=[0.75	1.15	2	2.1	0.6	200	5.5	1.1	0.2	0.3	0.95	0.03	0.34	0.13	1.15 1.5	1100	1.2	35	5.5 1	10000	75	5.5];  %pessimistic
    
%% Variable initialization
    ming=zeros(length(titles),5);
    maxg=zeros(length(titles),5);
    tiledlayout(5,5);
    counter=1;
    highval=zeros(21,1);
    lowval=zeros(21,1);
    bc1=bc;
    lb1=lb;
    ub1=ub;
        
%% Single Variable Sensitivity Analysis
    for main=1:length(titles)
        
    %% determination of upper and lower bound values    
        lowbound=zeros(1,length(titles));
        upbound=zeros(1,length(titles));
        basecase=bc1;
        for op=1:length(titles)
            if op==main
                lowbound(op)=basecase(op)*(1+lb1(op));
                upbound(op)=basecase(op)*(1+ub1(op));
            else
                lowbound(op)=basecase(op);
                upbound(op)=basecase(op);
            end
        end

    %% parameter initialization (within SVSA loop)
        cy=2;                 % construction years
        cd=[0.6 0.4 0 0 0];   % CAPEX split between the construction years 
        
        pl=20;                % Plant life
        SF = 0.9;             % Stream factor
        tpa=2000000;          % desired C2H4 in MMT/yr or TPA
        tr=0.39;              % tax rate
        lc=10;                % Land cost $10 million 
        c2cconv=100;          % CO3 to CO2 conversion
        h2ocost=0.1;          % $/m3
        kohcost=0.0005;       % $/kg
    
    %% variable definition (within SVSA loop)
        Result=zeros(999,5);              % storing the results of all the sensitivity runs   
        corr=zeros(1000,length(titles));
        proval=zeros(1000,1);
        op=zeros(1000,1);
        ele=zeros(1000,1);
        fci=zeros(1000,1);
        soln=zeros(1000,1);
    
        capdis=cd; 
        depp = [0.1 0.18 0.14 0.12 0.09 0.07370 0.06550 0.06550 0.06550 0.06550 0.03]; % 10-yr MACRS: Modified Accelerated Cost Recovery System
        projLife = pl;
        cons_yrs = cy;
        desired_C2H4=tpa;
        
        Taxrate = tr; 
        landcost = lc;
        CWT = 0 ;
        WACC=0.1;           % Weight Averaged Cost of Capital
        Salval=0;
          
        totyrs = cons_yrs + projLife + 1;
        
        % Calculating Cummulative Prpbability Dist based on triangular dist formulas when x= b
        probDist=zeros(1,length(titles));
        for i = 1:length(titles)
            if lowbound(i) ~= upbound(i) && lowbound(i) ~= basecase(i)                   %avoiding 0 in denominator
                probDist(i) = (basecase(i) - lowbound(i)) / (upbound(i) - lowbound(i));
            end
        end
    
        num=bc1; %% num is the basecase
                
        
        %% CCU calculations
        capturerate = 31.501*num(1) -0.1837;                    % t-CO2/m2/yr : correlation from CCU_calc sheet in data folder
        kohconc = 2.5;                                          % mol/L
        RR = 1/num(9);                                          % approximate recycle ratio
        co22bcap= (desired_C2H4* 1000/ 28/ 365/ 86400) /num(1) /num(5) /num(9) /(num(13)/(num(14)+num(13))) *2 /SF; % kmol/s : desired_C2H4 / capt_eff / rejection / SP_conv / Selecitivity_C2H4 (!= FE_C2H4) *2 carbonates per mol of C2+ prod
        
        if num(5)<0.75
            co22bcap = co22bcap /num(5)*RR;                       % lower carbonate rejection (<75%) requires 2-stage upstream membrane concentrator unit (pessimistic case)
        end
        
        kohCCU= co22bcap/0.13*56*RR;                            % kg/s, assuming 0.13 mol CO2/mol KOH uptake (Holmes and Keith, Phy. Trans. R.S.C 2012)
        waterCCU= kohCCU/56*1000/kohconc+2*18*co22bcap;         % kg/s
        CCUarea= co22bcap*44/capturerate/1000*3600*24*365*SF;   % m^2
        FCILCCU= 4668*(CCUarea^0.85)*num(2)*1.18*1.5;
        
        CRMCCU = (((waterCCU)*h2ocost*0.001)+(kohCCU*kohcost))*365*24*3600*SF/1000000;
        CRRCCU = FCILCCU*0.05/1000000;
        CUTCCU = (SF*365*24)*125.44*1.6/(60/100)/1000*CCUarea*num(12)/1000000; % packing pr drop = 125.44 Pa/m, air velocity =1.6 m/s, fan efficiency = 60%
        
        COMdCCU = 0.18 * (FCILCCU/1000000) + 1.23* (CUTCCU + CRMCCU + CRRCCU);
        
        CRF = 12.5/100;
        lvd_cost = (FCILCCU*CRF+COMdCCU*1000000)/(capturerate*CCUarea)/SF;

        %% GO Membrane concentrator calculations
        % CO3 concentrating stages
        if num(5)>0.8
            rej_co3= num(1)* co22bcap* num(5)* RR;              % kmol/s
            waterperm= waterCCU -  rej_co3*(1000)/num(3);       % kg/s (assuming density =1)
            TMP = 50;                                           % bar : transmembrane pressure
            memarea_t= waterperm/num(4)/TMP*3600;               % m^2 
        else
            rej_co3_1= num(1)* co22bcap* num(5)*RR;             % kmol/s
            waterperm_1= waterCCU*0.5 -  rej_co3_1*(1000)/1.25; % kg/s
            TMP = 50;                                           % bar : transmembrane pressure
            memarea_1= waterperm_1/num(4)/TMP*3600;             % m^2 
        
            rej_co3_2= num(1)* co22bcap* num(5) *0.45 *RR;       % kmol/s
            waterperm_2= waterCCU*0.5 -  rej_co3_2*(1000)/num(3);% kg/s
            TMP = 50;                                            % bar : transmembrane pressure
            rej_co3 = rej_co3_1 + rej_co3_2;
            memarea_t= waterperm_2/num(4)/TMP*3600 + memarea_1;  % m^2
            waterperm = (waterperm_1 + waterperm_2)*0.25;
        end
        
        % Recovery stage
        rejco3_r= (co22bcap*RR - rej_co3)* num(5);               % kmol/s
        waterperm_r= waterperm -  rejco3_r*(1000)/0.5;           % kg/s
        TMP = 50;                                                % bar : transmembrane pressure
        memarea_r= waterperm_r/num(4)/TMP*3600;                  % m^2 
        
        FCILGO=num(6)*((memarea_t + memarea_r)^0.85)*(pl/num(7))*1.5*1.18; 
        CUTGO=(waterCCU)*3.6*1000* 9.81* 1*SF*24*365*num(12)/3.6e6/1000000; % All in SI: Q* rho* g* head* eleccost/ 3.6E6
        CRRGO = FCILGO*0.05/1000000;
        
        %% electrolyzer calc
        n_C2H4=12;
        FECO  =0;
        FECH4 =0;
        FEH2  =1- num(13) -num(14) -FECO -FECH4;
        V = 4.26;         % Calculated using the simple voltammetric model described in the SI

        max_co2 = co22bcap *num(1) *num(5) *RR;      % kmol/s : max CO2 to the EchemR to convert to product
        Itot = (desired_C2H4* 1000/ 28/ 365/ 86400)* n_C2H4 *96485 *1000 /num(13)/SF ; % A
        Elecarea = Itot/num(10)/10000;               % m^2
        power = Itot * V/1000;                       % in kW 
        
        % stack design parameters
        cellperst = 150;
        areapercel = num(16);                        % m^2
        areaperst = areapercel* cellperst;           % m^2
        stacktot = round(Elecarea/areaperst);      
        
        % stack costing: refer stack costing excel sheet
        part_coeff = [194.7371361, 95.1714818, 83.5651965, 135.5265756, 257.659476, 29.3912632, 239.0894003];   % BPM w/ catholyte 100 from Badgett et al. J. of cleaner prod
        
        cath_flow = part_coeff(1)*(areaperst/0.0961) ^0.85/1000000*stacktot;
        cath_gdl  = part_coeff(2)*(areaperst/0.0676) ^0.85/1000000*stacktot;
        bpm       = part_coeff(3)*(areaperst/0.0961) ^0.85/1000000*stacktot;
        anod_gdl  = part_coeff(4)*(areaperst/0.0676) ^0.85/1000000*stacktot;
        mea_fr    = part_coeff(5)*(areaperst/0.0961) ^0.85/1000000*stacktot;
        bpp       = part_coeff(6)*(areaperst/0.0961) ^0.85/1000000*stacktot;
        stack_ass = part_coeff(7)*(areaperst/0.0961) ^0.85/1000000*stacktot;
        FCILEchemR  = (cath_flow + cath_gdl + bpm + anod_gdl + mea_fr + bpp + stack_ass)*1000000*num(8)*1.5*1.18;
        CUTEchemR = power*num(12)*24*365*SF/1000000;
        
        % Product distribution 
        h2   = Itot*FEH2*3600/1000*2.015/96485/2*SF;        % kg/hr
        CO   = Itot*FECO*3600/1000*28.01/96485/2*SF;        % kg/hr
        ethy = Itot*num(13)*3600/1000*28/96485/12*SF;       % kg/hr
        etha = Itot*num(14)*3600/1000*46/96485/12*SF;       % kg/hr 
        O2   = Itot*3600/1000/4/96485*32;                % kg/hr 
        H2Oin  = (ethy/28.05*4+etha/46.06*3+CO/28.01)*18 + (rej_co3*(1000)/num(3)*3600);               % kg/hr 
        H2Oreq = ((ethy/28*12)+(h2/2*2)+(etha/46*12)+(CO/28*2)+(ethy/28*2 + etha/46*2 + CO/28*1))*18;  % kg/hr 
        
        CRMEchemR = (H2Oreq)*h2ocost/1000*365*24*SF/1000000;
        CRREchemR = FCILEchemR* 0.05 /1000000;
        
        %% Downstream separations
        % cathode gas separation
        x_H2 = h2/2/(h2/2+ethy/28);                                      % feed H2 mol fraction
        Pf = 2;                                                          % bar : feed pressure
        ZHFarea=(h2/2/3600)/(num(17)*3.35e-10)/(Pf* 100000* x_H2)*1000;  % values subject to changes
        FCILZHF=num(19)*((ZHFarea)^0.85)*(pl/num(20))*1.5*1.18; 
        CUTZHF =(h2)*1* 9.81* 2*num(12)*24*365/1000000; % H2 density = 0.08375 kg/m^3
        CRRZHF = FCILZHF* 0.05 /1000000;
        
        % cathode liq separation only for optimistic case. The actual costs
        % for the column is obtained from Aspen Economic analyzer and is
        % part of the BoP costs. Added here to account for the separation
        % cost estimation.
        
        % (1) distillation
        flow_in = (H2Oin-H2Oreq + etha)/60;                             % L/min, density =1 
        FCILdis = 4162240*((flow_in/1000)^0.7)*1.5*1.18;                % cost calc from Jouny et al.    
        CUTdis=(9895*(flow_in/1000)^0.7 * 365*SF )/1000000;             % coolant flow+ heating requirement
        
        % (2) pervaporation
        q_h2o = etha*0.05;                                             % water in feed
        Pf = 1;                                                        % bar : feed pressure
        pervaparea= (q_h2o*0.999)/(num(21)/18);                        % desired ethanol purity in permeate
        FCILpervap= num(23)*((pervaparea)^0.85)*(pl/num(24))*1.5*1.18; 
        CUTpervap = q_h2o*1* 9.81* 1*num(12)*24*365/1000000; 
        CRRpervap = FCILpervap* 0.05 /1000000;
        
        
        %% total CAPEX & OPEX
        CUTI = CUTCCU + CUTGO + CUTEchemR + CUTZHF + CUTdis + CUTpervap;      % million USD
        FCILbop = 0.35*(FCILEchemR/0.65);                                     % USD 
        FCIL = FCILCCU + FCILGO + FCILEchemR + FCILZHF + FCILdis + FCILbop + FCILpervap;   % USD     
        CRR = CRRCCU + CRRGO + CRREchemR + CRRZHF + CRRpervap;                % million USD: Running and replacement
        CRM = CRMCCU + CRMEchemR;                                             % million USD: Raw materials
        NOL = (6.29 + 31.7* (0)^2 + 0.23*(8+stacktot))^0.5;                   % Eqn 8.3 richard turton, 8 unit ops
        COL = 4.5* NOL* (32.77*2000)/1000000;                                 % million USD: 245 shifts per operator, 1095 operating shifts per plant: 4.5 operators per required operator for the plant; per operator cost is 32.77 usd/hr, 2000 hrs per yr: https://www.bls.gov/regions/southwest/news-release/employercostsforemployeecompensation_regions.htm 
        workex = 0.1*CRM + 0.1*FCIL/1000000 + 0.1*COL;                        % working capital (millions USD)
        COMd = 0.18 * (FCIL/1000000) + 2.76 * COL + 1.23* (CUTI + CWT + CRM + CRR); % million USD
        
        %% Ethylene base case production cost (bcpc) estimation
        p2=zeros(20000,1);  
        j=1;
        for bcpc=0:0.001:200                                     
            prodincome=((0.6*CO + 3.9*h2 + 0.1345*O2 + 1*etha)*num(11) + bcpc*ethy )*24*365 /1000000;
        
            %% Profitability analysis
            % defining arrays
            Arr=zeros(100,5);
        
            %year 0 done separately as b4
            Start = 1;
            
            Arr(Start, 1) = -landcost;              % Non-discounted After tax cash flow
            Arr(Start, 2) = Arr(Start, 1);          % Cummulative Non-discounted After tax cash flow
            Arr(Start, 3) = Arr(Start, 1);          % Cummulative Discounted After tax cash flow
        
            Start = 2;
            
            while Start <= cons_yrs+1
                Arr(Start, 1) = -FCIL / 10 ^ 6 * capdis(Start - 1);
                if Start == cons_yrs+1 
                    Arr(Start, 1) = -(workex) - ((FCIL / 10 ^ 6) * capdis(Start - 1));
                end
                Arr(Start, 2) = Arr(Start, 1) + Arr(Start - 1, 2);
                Arr(Start, 3) = Arr(Start, 1) / ((1 + WACC) ^ (Start-1)) + Arr(Start - 1, 3);
                Start = Start + 1;
            end
            t=Start;        
            while t <= totyrs
                if (Start-cons_yrs>11)
                    depp(Start-cons_yrs)=0;
                end
                Arr(Start, 1) = (prodincome - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs-1)) * (1 - Taxrate) + (FCIL / 10 ^ 6) * depp(Start - cons_yrs-1);
                Arr(Start, 2) = Arr(Start, 1) + Arr(Start - 1, 2);
                Arr(Start, 3) = Arr(Start, 1) / ((WACC + 1) ^ (Start-1)) + Arr(Start - 1, 3);
                Arr(Start, 4) = (prodincome - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs)) * (1 - Taxrate);
                
                if Start == totyrs 
                    Arr(Start, 1) = (prodincome + Salval / 10 ^ 6 - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs)) * (1 - Taxrate) + (FCIL / 10 ^ 6) * depp(Start - cons_yrs) + workex + landcost;
                    Arr(Start, 2) = Arr(Start, 1) + Arr(Start - 1, 2);
                    Arr(Start, 3) = Arr(Start, 1) / ((WACC + 1) ^ (Start-1)) + Arr(Start - 1, 3);
                    Arr(Start, 4) = (prodincome - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs)) * (1 - Taxrate);
                end
                t = t + 1;
                Start = Start + 1;
            end
            
        %% (1) NPV value calculation
            Start = totyrs;
            NPV = Arr(Start, 3);
        
            if j~=1 && NPV>0  
                break
            end
            j=j+1;
        end
        
    % montecarlo simulation()
         for i = 1:1000
            %Defining random numbers for MC simulation
            randa=rand(1,length(titles));
            num=zeros(1,length(titles));
        
            %Calculation NPV with given Probability Dist(rando value above) based on triangular dist formulas
            for t = 1:length(titles)
                if upbound(t) == lowbound(t)
                    num(t) = upbound(t);
                elseif randa(t) <= probDist(t)
                    num(t) = (randa(t) * (upbound(t) - lowbound(t)) * (basecase(t) - lowbound(t))) ^ 0.5 + lowbound(t);   %using x<b formula
                else
                    %using x>b fomula, which gives a quad eqn: solving using the quadratic formula
                    %ax^2+bX+c=0
                    c = ((randa(t)-(basecase(t) - lowbound(t)) / (upbound(t) - lowbound(t))) * (upbound(t) - lowbound(t)) * (-basecase(t) + upbound(t))) - basecase(t) ^ 2 + 2 * basecase(t) * upbound(t);
                    b = 2 * upbound(t);
                    
                    root1 = (b + (b ^ 2 - 4 * c) ^ 0.5) / 2;
                    root2 = (b - (b ^ 2 - 4 * c) ^ 0.5) / 2;
        
                    if root1 > upbound(t) 
                        num(t) = root2;
                    else 
                        num(t) = root1;
                    end
                end
            end
               
            % Same calc as above, avoided defn as a function
            %% CCU calculations
            capturerate = 31.501*num(1) -0.1837;                    % t-CO2/m2/yr : correlation from CCU_calc sheet in data folder
            kohconc = 2.5;                                          % mol/L
            RR = 1/num(9);                                               % recycle ratio
            co22bcap= (desired_C2H4* 1000/ 28/ 365/ 86400) /num(1) /num(5) /num(9) /(num(13)/(num(14)+num(13))) *2 /RR; % kmol/s : desired_C2H4 / capt_eff / rejection / SP_conv / Selecitivity_C2H4 (!= FE_C2H4) *2 carbonates per mol of C2+ prod
            
            if num(5)<0.75
                co22bcap = co22bcap /0.45*RR;                          % 45% rej for stage 2
            end
            
            kohCCU= co22bcap/0.13*56*RR;                            % kg/s, assuming 0.13 mol CO2/mol KOH uptake
            waterCCU= kohCCU/56*1000/kohconc;                       % kg/s
            CCUarea= co22bcap*44/capturerate/1000*3600*24*365*SF;   % m^2
            
            FCILCCU= 4668*(CCUarea^0.85)*num(2)*1.18*1.5;
            CRMCCU = (((waterCCU)*h2ocost*0.001)+(kohCCU*kohcost))*365*24*3600*SF/1000000;
            CRRCCU = FCILCCU*0.05/1000000;
            CUTCCU = (SF*365*24)*125.44*1.6/(60/100)/1000*CCUarea*num(12)/1000000; % packing pr drop = 125.44 Pa/m, air velocity =1.6 m/s, fan efficiency = 60%
            
            COMdCCU = 0.18 * (FCILCCU/1000000) + 1.23* (CUTCCU + CRMCCU + CRRCCU);
            
            CRF = 12.5/100;
            lvd_cost = (FCILCCU*CRF+COMdCCU*1000000)/(capturerate*CCUarea)/SF;

            %% GO Membrane concentrator calculations
            % CO3 concentrating stages
            if num(5)>0.8
                rej_co3= num(1)* co22bcap* num(5)* RR;              % kmol/s
                waterperm= waterCCU -  rej_co3*(1000)/num(3);       % kg/s (assuming density =1)
                TMP = 50;                                           % bar : transmembrane pressure
                memarea_t= waterperm/num(4)/TMP*3600;               % m^2 
            else
                rej_co3_1= num(1)* co22bcap* num(5)*RR;             % kmol/s
                waterperm_1= waterCCU*0.5 -  rej_co3_1*(1000)/1.25; % kg/s
                TMP = 50;                                           % bar : transmembrane pressure
                memarea_1= waterperm_1/num(4)/TMP*3600;             % m^2 
            
                rej_co3_2= num(1)* co22bcap* num(5) *0.45 *RR;       % kmol/s
                waterperm_2= waterCCU*0.5 -  rej_co3_2*(1000)/num(3);% kg/s
                TMP = 50;                                            % bar : transmembrane pressure
                rej_co3 = rej_co3_1 + rej_co3_2;
                memarea_t= waterperm_2/num(4)/TMP*3600 + memarea_1;  % m^2
                waterperm = (waterperm_1 + waterperm_2)*0.25;
            end
            
            % Recovery stage
            rejco3_r= (co22bcap*RR - rej_co3)* num(5);               % kmol/s
            waterperm_r= waterperm -  rejco3_r*(1000)/0.5;           % kg/s
            TMP = 50;                                                % bar : transmembrane pressure
            memarea_r= waterperm_r/num(4)/TMP*3600;                  % m^2 
            
            FCILGO=num(6)*((memarea_t + memarea_r)^0.85)*(pl/num(7))*1.5*1.18; 
            CUTGO=(waterCCU)*3.6*1000* 9.81* 1*SF*24*365*num(12)/3.6e6/1000000; % All in SI: Q* rho* g* head* eleccost/ 3.6E6
            CRRGO = FCILGO*0.05/1000000;
            
            %% electrolyzer calc
            n_C2H4=12;
            FECO  =0;
            FECH4 =0;
            FEH2  =1- num(13) -num(14) -FECO -FECH4;
            V = (4.1* num(10)+ 2.3943)* num(15);         % Full cell voltage in volts from shin et al. (BPM), operating in the linear range of the Tafel plot
            % V = (0.2508* ln(num(10))+ 1.106)* num(15)  % Full cell voltage in volts from shin et al. (AEM)
            max_co2 = co22bcap *num(1) *num(5) *RR;      % kmol/s : max CO2 to the EchemR to convert to product
            Itot = (desired_C2H4* 1000/ 28/ 365/ 86400)* n_C2H4 *96485 *1000 /num(13)/SF ; % A
            Elecarea = Itot/num(10)/10000;               % m^2
            power = Itot * V/1000;                       % in kW 
            
            % stack design parameters
            cellperst = 150;
            areapercel = num(16);                        % m^2
            areaperst = areapercel* cellperst;           % m^2
            stacktot = round(Elecarea/areaperst);      
            
            % stack costing: refer stack costing excel sheet
            part_coeff = [194.7371361, 95.1714818, 83.5651965, 135.5265756, 257.659476, 29.3912632, 239.0894003];   % BPM w/ catholyte 100
            
            cath_flow = part_coeff(1)*(areaperst/0.0961) ^0.85/1000000*stacktot;
            cath_gdl  = part_coeff(2)*(areaperst/0.0676) ^0.85/1000000*stacktot;
            bpm       = part_coeff(3)*(areaperst/0.0961) ^0.85/1000000*stacktot;
            anod_gdl  = part_coeff(4)*(areaperst/0.0676) ^0.85/1000000*stacktot;
            mea_fr    = part_coeff(5)*(areaperst/0.0961) ^0.85/1000000*stacktot;
            bpp       = part_coeff(6)*(areaperst/0.0961) ^0.85/1000000*stacktot;
            stack_ass = part_coeff(7)*(areaperst/0.0961) ^0.85/1000000*stacktot;
            FCILEchemR  = (cath_flow + cath_gdl + bpm + anod_gdl + mea_fr + bpp + stack_ass)*1000000*num(8)*1.5*1.18;
            CUTEchemR = power*num(12)*24*365*SF/1000000;
            
            % Product distribution 
            h2   = Itot*FEH2*3600/1000*2.015/96485/2*SF;        % kg/hr
            CO   = Itot*FECO*3600/1000*28.01/96485/2*SF;        % kg/hr
            ethy = Itot*num(13)*3600/1000*28/96485/12*SF;       % kg/hr
            etha = Itot*num(14)*3600/1000*46/96485/12*SF;       % kg/hr 
            O2   = Itot*3600/1000/4/96485*32;                % kg/hr 
            H2Oin  = (ethy/28.05*4+etha/46.06*3+CO/28.01)*18 + (rej_co3*(1000)/num(3)*3600);               % kg/hr 
            H2Oreq = ((ethy/28*12)+(h2/2*2)+(etha/46*12)+(CO/28*2)+(ethy/28*2 + etha/46*2 + CO/28*1))*18;  % kg/hr 
            
            CRMEchemR = (H2Oreq)*h2ocost/1000*365*24*SF/1000000;
            CRREchemR = FCILEchemR* 0.05 /1000000;
            
            %% Downstream separations
            % cathode gas separation
            x_H2 = h2/2/(h2/2+ethy/28);                                      % feed H2 mol fraction
            Pf = 2;                                                          % bar : feed pressure
            ZHFarea=(h2/2/3600)/(num(17)*3.35e-10)/(Pf* 100000* x_H2)*1000;  % values subject to changes
            FCILZHF=num(19)*((ZHFarea)^0.85)*(pl/num(20))*1.5*1.18; 
            CUTZHF =(h2)*1* 9.81* 2*num(12)*24*365/1000000; % H2 density = 0.08375 kg/m^3
            CRRZHF = FCILZHF* 0.05 /1000000;
            
            % cathode liq separation
            % (1) distillation
            flow_in = (H2Oin-H2Oreq + etha)/60;                 % L/min, density =1 
            FCILdis = 4162240*((flow_in/1000)^0.7)*1.5*1.18;    % cost calc from Jouny et al.    
            CUTdis=(9895*(flow_in/1000)^0.7 * 365*SF )/1000000;       % coolant flow+ heating requirement
            
            % (2) pervaporation
            q_h2o = etha*0.05;                                             % water in feed
            Pf = 1;                                                        % bar : feed pressure
            pervaparea= (q_h2o*0.999)/(num(21)/18);                        % desired ethanol purity in permeate
            FCILpervap= num(23)*((pervaparea)^0.85)*(pl/num(24))*1.5*1.18; 
            CUTpervap = q_h2o*1* 9.81* 1*num(12)*24*365/1000000; 
            CRRpervap = FCILpervap* 0.05 /1000000;
            
            
            %% total CAPEX & OPEX
            CUTI = CUTCCU + CUTGO + CUTEchemR + CUTZHF + CUTdis + CUTpervap;      % million USD
            FCILbop = 0.35*(FCILEchemR/0.65);                                     % USD 
            FCIL = FCILCCU + FCILGO + FCILEchemR + FCILZHF + FCILdis + FCILbop + FCILpervap;   % USD     
            CRR = CRRCCU + CRRGO + CRREchemR + CRRZHF + CRRpervap;                % million USD: Running and replacement
            CRM = CRMCCU + CRMEchemR;                                             % million USD: Raw materials
            NOL = (6.29 + 31.7* (0)^2 + 0.23*(8+stacktot))^0.5;                   % Eqn 8.3 richard turton, 8 unit ops
            COL = 4.5* NOL* (32.77*2000)/1000000;                                 % million USD: 245 shifts per operator, 1095 operating shifts per plant: 4.5 operators per required operator for the plant; per operator cost is 32.77 usd/hr, 2000 hrs per yr: https://www.bls.gov/regions/southwest/news-release/employercostsforemployeecompensation_regions.htm 
            workex = 0.1*CRM + 0.1*FCIL/1000000 + 0.1*COL;                        % working capital (millions USD)
            COMd = 0.18 * (FCIL/1000000) + 2.76 * COL + 1.23* (CUTI + CWT + CRM + CRR); % million USD
            
            %% Ethylene base case production cost (bcpc) estimation
            p2=zeros(20000,1);  
            j=1;
            for msp=0:0.002:200                                     
                prodincome=((0.6*CO + 3.9*h2 + 0.1345*O2 + 1*etha)*num(11) + msp*ethy )*24*365 /1000000;
            
                %% Profitability analysis
                % defining arrays
                Arr=zeros(100,5);
            
                %year 0 done separately as b4
                Start = 1;
                
                Arr(Start, 1) = -landcost;              % Non-discounted After tax cash flow
                Arr(Start, 2) = Arr(Start, 1);          % Cummulative Non-discounted After tax cash flow
                Arr(Start, 3) = Arr(Start, 1);          % Cummulative Discounted After tax cash flow
            
                Start = 2;
                
                while Start <= cons_yrs+1
                    Arr(Start, 1) = -FCIL / 10 ^ 6 * capdis(Start - 1);
                    if Start == cons_yrs+1 
                        Arr(Start, 1) = -(workex) - ((FCIL / 10 ^ 6) * capdis(Start - 1));
                    end
                    Arr(Start, 2) = Arr(Start, 1) + Arr(Start - 1, 2);
                    Arr(Start, 3) = Arr(Start, 1) / ((1 + WACC) ^ (Start-1)) + Arr(Start - 1, 3);
                    Start = Start + 1;
                end
                t=Start;        
                while t <= totyrs
                    if (Start-cons_yrs>11)
                        depp(Start-cons_yrs)=0;
                    end
                    Arr(Start, 1) = (prodincome - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs-1)) * (1 - Taxrate) + (FCIL / 10 ^ 6) * depp(Start - cons_yrs-1);
                    Arr(Start, 2) = Arr(Start, 1) + Arr(Start - 1, 2);
                    Arr(Start, 3) = Arr(Start, 1) / ((WACC + 1) ^ (Start-1)) + Arr(Start - 1, 3);
                    Arr(Start, 4) = (prodincome - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs)) * (1 - Taxrate);
                    
                    if Start == totyrs 
                        Arr(Start, 1) = (prodincome + Salval / 10 ^ 6 - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs)) * (1 - Taxrate) + (FCIL / 10 ^ 6) * depp(Start - cons_yrs) + workex + landcost;
                        Arr(Start, 2) = Arr(Start, 1) + Arr(Start - 1, 2);
                        Arr(Start, 3) = Arr(Start, 1) / ((WACC + 1) ^ (Start-1)) + Arr(Start - 1, 3);
                        Arr(Start, 4) = (prodincome - COMd - (FCIL / 10 ^ 6) * depp(Start - cons_yrs)) * (1 - Taxrate);
                    end
                    t = t + 1;
                    Start = Start + 1;
                end
                
            %% (1) NPV value calculation
                Start = totyrs;
                NPV = Arr(Start, 3);
            
                if j~=1 && NPV>0  
                    break
                end
                j=j+1;
            end
            soln(i)=msp;
            corr(i,main)=num(main);
            proval(i)=prodincome;
            op(i)=COMd;
            ele(i)=FCILCCU;
            fci(i)=FCIL;
            
            %% (2) CCP value caclulation (in million USD)
            Result(i, 2) = Arr(Start, 2);
            
            %% (3) ROROI value calculation
            Result(i, 3) = mean(Arr(:,4)) / (FCIL / 10 ^ 6) *100;
            
            %% (4) DPBP value calculation
            Start = cons_yrs+1;
            capdisr = Start;
            
            FCILatStart = 0;
            
            for b = 1:cons_yrs+1
                FCILatStart = (FCILatStart + FCIL * capdis(b)) / (1 + WACC) ^ (b-1);
            end
                
            while Arr(Start, 3) <= (Arr(capdisr, 2) + FCILatStart / 10 ^ 6) && Arr(Start, 3) ~= 0
                Start = Start + 1;
            end
            
            %nth year plus the decimal amount required in the year by just taking percent FCI acquired
            if Start > totyrs 
                Result(i, 4) = 0; 
            else
                Result(i, 4) = Start - 1 - capdisr + ((Arr(Start - 1, 3) - (Arr(capdisr, 3) + FCIL / 10 ^ 6)) /(Arr(Start - 1, 3) - Arr(Start, 3)));
            end
            
            %% PBP value calculation
            Start = cons_yrs+1;
            capdisr = Start;
                
            while Arr(Start, 2) <= (Arr(capdisr, 2) + FCIL / 10 ^ 6) && Arr(Start, 1) ~= 0
                Start = Start + 1;
            end
            
            if Start > totyrs
                Result(i, 5) = Start; 
            else
                Result(i, 5) = Start - 1 - capdisr + (Arr(Start, 1) - (Arr(capdisr, 1) + FCIL / 10 ^ 6)) / (Arr(Start, 1) - Arr(Start-1, 1));
            end
        end
        %% writing to excel
        %SVSA
        low = min(soln);
        high = max(soln);
        pcval=zeros(21,1);
        pcBin=zeros(21,1);
        for l2 = 1:21
            pcval(l2)= low + (high - low) / 20 * (l2-1);
            for l1 = 1:1000
                if soln(l1) >= (low + (high - low) / 20 * (l2-1)) && soln(l1) <= (low + (high - low) / 20 * (l2))
                    pcBin(l2) = pcBin(l2) + 1;
                end
            end
        end
        pcvalc(:,counter)=pcval;
        pcvalc(:,counter+1)=pcBin/1000;
        counter=counter+2;
        nexttile
        y=pcBin/1000;
        ming(main)=min(pcval);
        maxg(main)=max(pcval);
        if length(unique(y))>1
            b = bar(pcval,y, 'FaceColor','flat');
            b.FaceColor = '#A2142F';
            title(titles(main))
            xlabel('MSP C2H4(USD/kg)')
            ylabel('p(x)')
        else 
            b = bar(msp,1,0.05, 'FaceColor','flat');
            b.FaceColor = '#A2142F';
            title(titles(main))
            xlabel('MSP C2H4(USD/kg)')
            ylabel('p(x)')
        end
    end    
    
    %% writing to excel
    %profitability
    l=zeros(5,1);
    h=zeros(5,1);
    iBin=zeros(5,length(titles));
    upval=zeros(5,length(titles));
    for g = 1:5
        low = min(Result(:, g));
        high = max(Result(:, g));
        l(g)=low;
        h(g)=high;
     
        for l2 = 1:11
            upval(g,l2) = low + (high - low) / 10 * (l2-1);
            for l1 = 1:1000
                if Result(l1, g) >= (low + (high - low) / 10 * (l2-1)) && Result(l1, g) <= (low + (high - low) / 10 * (l2))
                    iBin(g,l2) = iBin(g,l2) + 1;
                end
            end
        end
    end
    ff=cat(2,ming.',maxg.',l,h,iBin,upval);
FCILarr = [FCILCCU, FCILGO, FCILEchemR, FCILZHF, FCILdis, FCILpervap, FCILbop]/1000000;
% end
