  function targMap = targDataMap(),

  ;%***********************
  ;% Create Parameter Map *
  ;%***********************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 2;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc paramMap
    ;%
    paramMap.nSections           = nTotSects;
    paramMap.sectIdxOffset       = sectIdxOffset;
      paramMap.sections(nTotSects) = dumSection; %prealloc
    paramMap.nTotData            = -1;
    
    ;%
    ;% Auto data (Control_linear_P)
    ;%
      section.nData     = 75;
      section.data(75)  = dumData; %prealloc
      
	  ;% Control_linear_P.GT400SVInitialization1_P1
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_P.TransportDelay1_Delay
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Control_linear_P.TransportDelay1_InitOutput
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Control_linear_P.TransportDelay_Delay
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Control_linear_P.TransportDelay_InitOutput
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% Control_linear_P.radtodeg1_Gain
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% Control_linear_P.Saturation_UpperSat
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% Control_linear_P.Saturation_LowerSat
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 7;
	
	  ;% Control_linear_P.GetCurrentAxisPosition_P1
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 8;
	
	  ;% Control_linear_P.Gain_Gain
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 9;
	
	  ;% Control_linear_P.Constant7_Value
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 10;
	
	  ;% Control_linear_P.radtodeg3_Gain
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 11;
	
	  ;% Control_linear_P.Constant_Value
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 12;
	
	  ;% Control_linear_P.Step1_Time
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 13;
	
	  ;% Control_linear_P.Step1_Y0
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 14;
	
	  ;% Control_linear_P.Step1_YFinal
	  section.data(16).logicalSrcIdx = 15;
	  section.data(16).dtTransOffset = 15;
	
	  ;% Control_linear_P.TransferFcn_A
	  section.data(17).logicalSrcIdx = 16;
	  section.data(17).dtTransOffset = 16;
	
	  ;% Control_linear_P.TransferFcn_C
	  section.data(18).logicalSrcIdx = 17;
	  section.data(18).dtTransOffset = 18;
	
	  ;% Control_linear_P.TransportDelay1_Delay_m
	  section.data(19).logicalSrcIdx = 20;
	  section.data(19).dtTransOffset = 20;
	
	  ;% Control_linear_P.TransportDelay1_InitOutput_g
	  section.data(20).logicalSrcIdx = 21;
	  section.data(20).dtTransOffset = 21;
	
	  ;% Control_linear_P.TransportDelay_Delay_b
	  section.data(21).logicalSrcIdx = 22;
	  section.data(21).dtTransOffset = 22;
	
	  ;% Control_linear_P.TransportDelay_InitOutput_m
	  section.data(22).logicalSrcIdx = 23;
	  section.data(22).dtTransOffset = 23;
	
	  ;% Control_linear_P.radtodeg1_Gain_b
	  section.data(23).logicalSrcIdx = 24;
	  section.data(23).dtTransOffset = 24;
	
	  ;% Control_linear_P.Saturation_UpperSat_b
	  section.data(24).logicalSrcIdx = 25;
	  section.data(24).dtTransOffset = 25;
	
	  ;% Control_linear_P.Saturation_LowerSat_m
	  section.data(25).logicalSrcIdx = 26;
	  section.data(25).dtTransOffset = 26;
	
	  ;% Control_linear_P.GetCurrentAxisPosition2_P1
	  section.data(26).logicalSrcIdx = 27;
	  section.data(26).dtTransOffset = 27;
	
	  ;% Control_linear_P.Gain2_Gain
	  section.data(27).logicalSrcIdx = 28;
	  section.data(27).dtTransOffset = 28;
	
	  ;% Control_linear_P.radtodeg7_Gain
	  section.data(28).logicalSrcIdx = 29;
	  section.data(28).dtTransOffset = 29;
	
	  ;% Control_linear_P.TransportDelay2_Delay
	  section.data(29).logicalSrcIdx = 30;
	  section.data(29).dtTransOffset = 30;
	
	  ;% Control_linear_P.TransportDelay2_InitOutput
	  section.data(30).logicalSrcIdx = 31;
	  section.data(30).dtTransOffset = 31;
	
	  ;% Control_linear_P.Constant_Value_n
	  section.data(31).logicalSrcIdx = 32;
	  section.data(31).dtTransOffset = 32;
	
	  ;% Control_linear_P.Step1_Time_n
	  section.data(32).logicalSrcIdx = 33;
	  section.data(32).dtTransOffset = 33;
	
	  ;% Control_linear_P.Step1_Y0_i
	  section.data(33).logicalSrcIdx = 34;
	  section.data(33).dtTransOffset = 34;
	
	  ;% Control_linear_P.Step1_YFinal_f
	  section.data(34).logicalSrcIdx = 35;
	  section.data(34).dtTransOffset = 35;
	
	  ;% Control_linear_P.TransferFcn_A_b
	  section.data(35).logicalSrcIdx = 36;
	  section.data(35).dtTransOffset = 36;
	
	  ;% Control_linear_P.TransferFcn_C_g
	  section.data(36).logicalSrcIdx = 37;
	  section.data(36).dtTransOffset = 38;
	
	  ;% Control_linear_P.k_Gain
	  section.data(37).logicalSrcIdx = 40;
	  section.data(37).dtTransOffset = 40;
	
	  ;% Control_linear_P.radtodeg3_Gain_l
	  section.data(38).logicalSrcIdx = 41;
	  section.data(38).dtTransOffset = 42;
	
	  ;% Control_linear_P.GetCurrentAxisPosition1_P1
	  section.data(39).logicalSrcIdx = 42;
	  section.data(39).dtTransOffset = 43;
	
	  ;% Control_linear_P.Gain1_Gain
	  section.data(40).logicalSrcIdx = 43;
	  section.data(40).dtTransOffset = 44;
	
	  ;% Control_linear_P.radtodeg1_Gain_f
	  section.data(41).logicalSrcIdx = 44;
	  section.data(41).dtTransOffset = 45;
	
	  ;% Control_linear_P.TransferFcn_A_j
	  section.data(42).logicalSrcIdx = 45;
	  section.data(42).dtTransOffset = 46;
	
	  ;% Control_linear_P.TransferFcn_C_o
	  section.data(43).logicalSrcIdx = 46;
	  section.data(43).dtTransOffset = 48;
	
	  ;% Control_linear_P.TransferFcn_A_m
	  section.data(44).logicalSrcIdx = 49;
	  section.data(44).dtTransOffset = 50;
	
	  ;% Control_linear_P.TransferFcn_C_f
	  section.data(45).logicalSrcIdx = 50;
	  section.data(45).dtTransOffset = 52;
	
	  ;% Control_linear_P.radtodeg2_Gain
	  section.data(46).logicalSrcIdx = 53;
	  section.data(46).dtTransOffset = 54;
	
	  ;% Control_linear_P.radtodeg5_Gain
	  section.data(47).logicalSrcIdx = 54;
	  section.data(47).dtTransOffset = 62;
	
	  ;% Control_linear_P.radtodeg4_Gain
	  section.data(48).logicalSrcIdx = 55;
	  section.data(48).dtTransOffset = 63;
	
	  ;% Control_linear_P.Constant1_Value
	  section.data(49).logicalSrcIdx = 56;
	  section.data(49).dtTransOffset = 64;
	
	  ;% Control_linear_P.Constant2_Value
	  section.data(50).logicalSrcIdx = 57;
	  section.data(50).dtTransOffset = 65;
	
	  ;% Control_linear_P.Constant3_Value
	  section.data(51).logicalSrcIdx = 58;
	  section.data(51).dtTransOffset = 66;
	
	  ;% Control_linear_P.Constant4_Value
	  section.data(52).logicalSrcIdx = 59;
	  section.data(52).dtTransOffset = 67;
	
	  ;% Control_linear_P.Constant5_Value
	  section.data(53).logicalSrcIdx = 60;
	  section.data(53).dtTransOffset = 68;
	
	  ;% Control_linear_P.Gain_Gain_j
	  section.data(54).logicalSrcIdx = 61;
	  section.data(54).dtTransOffset = 69;
	
	  ;% Control_linear_P.Limit1_UpperSat
	  section.data(55).logicalSrcIdx = 62;
	  section.data(55).dtTransOffset = 70;
	
	  ;% Control_linear_P.Limit1_LowerSat
	  section.data(56).logicalSrcIdx = 63;
	  section.data(56).dtTransOffset = 71;
	
	  ;% Control_linear_P.Gain1_Gain_h
	  section.data(57).logicalSrcIdx = 64;
	  section.data(57).dtTransOffset = 72;
	
	  ;% Control_linear_P.Limit2_UpperSat
	  section.data(58).logicalSrcIdx = 65;
	  section.data(58).dtTransOffset = 73;
	
	  ;% Control_linear_P.Limit2_LowerSat
	  section.data(59).logicalSrcIdx = 66;
	  section.data(59).dtTransOffset = 74;
	
	  ;% Control_linear_P.radtodeg6_Gain
	  section.data(60).logicalSrcIdx = 67;
	  section.data(60).dtTransOffset = 75;
	
	  ;% Control_linear_P.Gain2_Gain_l
	  section.data(61).logicalSrcIdx = 68;
	  section.data(61).dtTransOffset = 76;
	
	  ;% Control_linear_P.Limit1_UpperSat_j
	  section.data(62).logicalSrcIdx = 69;
	  section.data(62).dtTransOffset = 77;
	
	  ;% Control_linear_P.Limit1_LowerSat_m
	  section.data(63).logicalSrcIdx = 70;
	  section.data(63).dtTransOffset = 78;
	
	  ;% Control_linear_P.Gain1_Gain_d
	  section.data(64).logicalSrcIdx = 71;
	  section.data(64).dtTransOffset = 79;
	
	  ;% Control_linear_P.Limit2_UpperSat_j
	  section.data(65).logicalSrcIdx = 72;
	  section.data(65).dtTransOffset = 80;
	
	  ;% Control_linear_P.Limit2_LowerSat_m
	  section.data(66).logicalSrcIdx = 73;
	  section.data(66).dtTransOffset = 81;
	
	  ;% Control_linear_P.TransferFcn2_A
	  section.data(67).logicalSrcIdx = 74;
	  section.data(67).dtTransOffset = 82;
	
	  ;% Control_linear_P.TransferFcn2_C
	  section.data(68).logicalSrcIdx = 75;
	  section.data(68).dtTransOffset = 83;
	
	  ;% Control_linear_P.SetCurrentAxisCommand1_P1
	  section.data(69).logicalSrcIdx = 78;
	  section.data(69).dtTransOffset = 84;
	
	  ;% Control_linear_P.TransferFcn1_A
	  section.data(70).logicalSrcIdx = 79;
	  section.data(70).dtTransOffset = 85;
	
	  ;% Control_linear_P.TransferFcn1_C
	  section.data(71).logicalSrcIdx = 80;
	  section.data(71).dtTransOffset = 86;
	
	  ;% Control_linear_P.SetCurrentAxisCommand2_P1
	  section.data(72).logicalSrcIdx = 83;
	  section.data(72).dtTransOffset = 87;
	
	  ;% Control_linear_P.Constant7_Value_l
	  section.data(73).logicalSrcIdx = 84;
	  section.data(73).dtTransOffset = 88;
	
	  ;% Control_linear_P.Constant1_Value_l
	  section.data(74).logicalSrcIdx = 85;
	  section.data(74).dtTransOffset = 89;
	
	  ;% Control_linear_P.Constant1_Value_d
	  section.data(75).logicalSrcIdx = 86;
	  section.data(75).dtTransOffset = 90;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
      section.nData     = 14;
      section.data(14)  = dumData; %prealloc
      
	  ;% Control_linear_P.ManualSwitch_CurrentSetting
	  section.data(1).logicalSrcIdx = 87;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_P.ManualSwitch5_CurrentSetting
	  section.data(2).logicalSrcIdx = 88;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Control_linear_P.ManualSwitch1_CurrentSetting
	  section.data(3).logicalSrcIdx = 89;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Control_linear_P.ManualSwitch4_CurrentSetting
	  section.data(4).logicalSrcIdx = 90;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Control_linear_P.ManualSwitch_CurrentSetting_e
	  section.data(5).logicalSrcIdx = 91;
	  section.data(5).dtTransOffset = 4;
	
	  ;% Control_linear_P.ManualSwitch5_CurrentSetting_d
	  section.data(6).logicalSrcIdx = 92;
	  section.data(6).dtTransOffset = 5;
	
	  ;% Control_linear_P.ManualSwitch1_CurrentSetting_a
	  section.data(7).logicalSrcIdx = 93;
	  section.data(7).dtTransOffset = 6;
	
	  ;% Control_linear_P.ManualSwitch4_CurrentSetting_k
	  section.data(8).logicalSrcIdx = 94;
	  section.data(8).dtTransOffset = 7;
	
	  ;% Control_linear_P.ManualSwitch2_CurrentSetting
	  section.data(9).logicalSrcIdx = 95;
	  section.data(9).dtTransOffset = 8;
	
	  ;% Control_linear_P.ManualSwitch1_CurrentSetting_f
	  section.data(10).logicalSrcIdx = 96;
	  section.data(10).dtTransOffset = 9;
	
	  ;% Control_linear_P.ManualSwitch2_CurrentSetting_l
	  section.data(11).logicalSrcIdx = 97;
	  section.data(11).dtTransOffset = 10;
	
	  ;% Control_linear_P.ManualSwitch3_CurrentSetting
	  section.data(12).logicalSrcIdx = 98;
	  section.data(12).dtTransOffset = 11;
	
	  ;% Control_linear_P.ManualSwitch2_CurrentSetting_m
	  section.data(13).logicalSrcIdx = 99;
	  section.data(13).dtTransOffset = 12;
	
	  ;% Control_linear_P.ManualSwitch3_CurrentSetting_d
	  section.data(14).logicalSrcIdx = 100;
	  section.data(14).dtTransOffset = 13;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(2) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (parameter)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    paramMap.nTotData = nTotData;
    


  ;%**************************
  ;% Create Block Output Map *
  ;%**************************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 1;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc sigMap
    ;%
    sigMap.nSections           = nTotSects;
    sigMap.sectIdxOffset       = sectIdxOffset;
      sigMap.sections(nTotSects) = dumSection; %prealloc
    sigMap.nTotData            = -1;
    
    ;%
    ;% Auto data (Control_linear_B)
    ;%
      section.nData     = 54;
      section.data(54)  = dumData; %prealloc
      
	  ;% Control_linear_B.Clock
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_B.FromWorkspace
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Control_linear_B.ManualSwitch5
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Control_linear_B.GetCurrentAxisPosition
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Control_linear_B.Gain
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% Control_linear_B.q1
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% Control_linear_B.radtodeg3
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% Control_linear_B.Sum2
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 7;
	
	  ;% Control_linear_B.FromWorkspace_b
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 8;
	
	  ;% Control_linear_B.ManualSwitch4
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 9;
	
	  ;% Control_linear_B.TransferFcn
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 10;
	
	  ;% Control_linear_B.Sum7
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 11;
	
	  ;% Control_linear_B.Clock_e
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 12;
	
	  ;% Control_linear_B.FromWorkspace_f
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 13;
	
	  ;% Control_linear_B.ManualSwitch5_k
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 14;
	
	  ;% Control_linear_B.GetCurrentAxisPosition2
	  section.data(16).logicalSrcIdx = 15;
	  section.data(16).dtTransOffset = 15;
	
	  ;% Control_linear_B.Gain2
	  section.data(17).logicalSrcIdx = 16;
	  section.data(17).dtTransOffset = 16;
	
	  ;% Control_linear_B.radtodeg7
	  section.data(18).logicalSrcIdx = 17;
	  section.data(18).dtTransOffset = 17;
	
	  ;% Control_linear_B.Sum1
	  section.data(19).logicalSrcIdx = 18;
	  section.data(19).dtTransOffset = 18;
	
	  ;% Control_linear_B.Fcn1
	  section.data(20).logicalSrcIdx = 19;
	  section.data(20).dtTransOffset = 19;
	
	  ;% Control_linear_B.FromWorkspace_bg
	  section.data(21).logicalSrcIdx = 20;
	  section.data(21).dtTransOffset = 20;
	
	  ;% Control_linear_B.ManualSwitch4_h
	  section.data(22).logicalSrcIdx = 21;
	  section.data(22).dtTransOffset = 21;
	
	  ;% Control_linear_B.TransferFcn_n
	  section.data(23).logicalSrcIdx = 22;
	  section.data(23).dtTransOffset = 22;
	
	  ;% Control_linear_B.Sum6
	  section.data(24).logicalSrcIdx = 23;
	  section.data(24).dtTransOffset = 23;
	
	  ;% Control_linear_B.radtodeg3_l
	  section.data(25).logicalSrcIdx = 24;
	  section.data(25).dtTransOffset = 24;
	
	  ;% Control_linear_B.GetCurrentAxisPosition1
	  section.data(26).logicalSrcIdx = 25;
	  section.data(26).dtTransOffset = 25;
	
	  ;% Control_linear_B.Gain1
	  section.data(27).logicalSrcIdx = 26;
	  section.data(27).dtTransOffset = 26;
	
	  ;% Control_linear_B.radtodeg1
	  section.data(28).logicalSrcIdx = 27;
	  section.data(28).dtTransOffset = 27;
	
	  ;% Control_linear_B.TransferFcn_k
	  section.data(29).logicalSrcIdx = 28;
	  section.data(29).dtTransOffset = 28;
	
	  ;% Control_linear_B.TransferFcn_l
	  section.data(30).logicalSrcIdx = 29;
	  section.data(30).dtTransOffset = 29;
	
	  ;% Control_linear_B.Sum8
	  section.data(31).logicalSrcIdx = 30;
	  section.data(31).dtTransOffset = 30;
	
	  ;% Control_linear_B.radtodeg2
	  section.data(32).logicalSrcIdx = 31;
	  section.data(32).dtTransOffset = 31;
	
	  ;% Control_linear_B.radtodeg5
	  section.data(33).logicalSrcIdx = 32;
	  section.data(33).dtTransOffset = 33;
	
	  ;% Control_linear_B.radtodeg4
	  section.data(34).logicalSrcIdx = 33;
	  section.data(34).dtTransOffset = 34;
	
	  ;% Control_linear_B.Sum9
	  section.data(35).logicalSrcIdx = 34;
	  section.data(35).dtTransOffset = 35;
	
	  ;% Control_linear_B.Sum10
	  section.data(36).logicalSrcIdx = 35;
	  section.data(36).dtTransOffset = 36;
	
	  ;% Control_linear_B.Limit1
	  section.data(37).logicalSrcIdx = 36;
	  section.data(37).dtTransOffset = 37;
	
	  ;% Control_linear_B.Limit2
	  section.data(38).logicalSrcIdx = 37;
	  section.data(38).dtTransOffset = 38;
	
	  ;% Control_linear_B.radtodeg6
	  section.data(39).logicalSrcIdx = 38;
	  section.data(39).dtTransOffset = 39;
	
	  ;% Control_linear_B.Limit1_a
	  section.data(40).logicalSrcIdx = 39;
	  section.data(40).dtTransOffset = 40;
	
	  ;% Control_linear_B.Limit2_o
	  section.data(41).logicalSrcIdx = 40;
	  section.data(41).dtTransOffset = 41;
	
	  ;% Control_linear_B.ManualSwitch2
	  section.data(42).logicalSrcIdx = 41;
	  section.data(42).dtTransOffset = 42;
	
	  ;% Control_linear_B.ManualSwitch1
	  section.data(43).logicalSrcIdx = 42;
	  section.data(43).dtTransOffset = 43;
	
	  ;% Control_linear_B.Sum1_p
	  section.data(44).logicalSrcIdx = 43;
	  section.data(44).dtTransOffset = 44;
	
	  ;% Control_linear_B.FromWorkspace_e
	  section.data(45).logicalSrcIdx = 44;
	  section.data(45).dtTransOffset = 45;
	
	  ;% Control_linear_B.ManualSwitch3
	  section.data(46).logicalSrcIdx = 45;
	  section.data(46).dtTransOffset = 46;
	
	  ;% Control_linear_B.Clock1
	  section.data(47).logicalSrcIdx = 46;
	  section.data(47).dtTransOffset = 47;
	
	  ;% Control_linear_B.FromWorkspace_c
	  section.data(48).logicalSrcIdx = 47;
	  section.data(48).dtTransOffset = 48;
	
	  ;% Control_linear_B.ManualSwitch3_o
	  section.data(49).logicalSrcIdx = 48;
	  section.data(49).dtTransOffset = 49;
	
	  ;% Control_linear_B.Clock1_i
	  section.data(50).logicalSrcIdx = 49;
	  section.data(50).dtTransOffset = 50;
	
	  ;% Control_linear_B.v1
	  section.data(51).logicalSrcIdx = 50;
	  section.data(51).dtTransOffset = 51;
	
	  ;% Control_linear_B.v2
	  section.data(52).logicalSrcIdx = 51;
	  section.data(52).dtTransOffset = 52;
	
	  ;% Control_linear_B.fs_d
	  section.data(53).logicalSrcIdx = 52;
	  section.data(53).dtTransOffset = 53;
	
	  ;% Control_linear_B.fd_d
	  section.data(54).logicalSrcIdx = 53;
	  section.data(54).dtTransOffset = 54;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(1) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (signal)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    sigMap.nTotData = nTotData;
    


  ;%*******************
  ;% Create DWork Map *
  ;%*******************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 6;
    sectIdxOffset = 1;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc dworkMap
    ;%
    dworkMap.nSections           = nTotSects;
    dworkMap.sectIdxOffset       = sectIdxOffset;
      dworkMap.sections(nTotSects) = dumSection; %prealloc
    dworkMap.nTotData            = -1;
    
    ;%
    ;% Auto data (Control_linear_DWork)
    ;%
      section.nData     = 5;
      section.data(5)  = dumData; %prealloc
      
	  ;% Control_linear_DWork.TransportDelay1_RWORK.modelTStart
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_DWork.TransportDelay_RWORK.modelTStart
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Control_linear_DWork.TransportDelay1_RWORK_l.modelTStart
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Control_linear_DWork.TransportDelay_RWORK_e.modelTStart
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Control_linear_DWork.TransportDelay2_RWORK.modelTStart
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 51;
      section.data(51)  = dumData; %prealloc
      
	  ;% Control_linear_DWork.TransportDelay1_PWORK.TUbufferPtrs
	  section.data(1).logicalSrcIdx = 5;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_DWork.TransportDelay_PWORK.TUbufferPtrs
	  section.data(2).logicalSrcIdx = 6;
	  section.data(2).dtTransOffset = 2;
	
	  ;% Control_linear_DWork.FromWorkspace_PWORK.TimePtr
	  section.data(3).logicalSrcIdx = 7;
	  section.data(3).dtTransOffset = 4;
	
	  ;% Control_linear_DWork.FromWorkspace_PWORK_c.TimePtr
	  section.data(4).logicalSrcIdx = 8;
	  section.data(4).dtTransOffset = 5;
	
	  ;% Control_linear_DWork.TransportDelay1_PWORK_i.TUbufferPtrs
	  section.data(5).logicalSrcIdx = 9;
	  section.data(5).dtTransOffset = 6;
	
	  ;% Control_linear_DWork.TransportDelay_PWORK_o.TUbufferPtrs
	  section.data(6).logicalSrcIdx = 10;
	  section.data(6).dtTransOffset = 8;
	
	  ;% Control_linear_DWork.FromWorkspace_PWORK_j.TimePtr
	  section.data(7).logicalSrcIdx = 11;
	  section.data(7).dtTransOffset = 10;
	
	  ;% Control_linear_DWork.TransportDelay2_PWORK.TUbufferPtrs
	  section.data(8).logicalSrcIdx = 12;
	  section.data(8).dtTransOffset = 11;
	
	  ;% Control_linear_DWork.FromWorkspace_PWORK_p.TimePtr
	  section.data(9).logicalSrcIdx = 13;
	  section.data(9).dtTransOffset = 13;
	
	  ;% Control_linear_DWork.Scope1_PWORK.LoggedData
	  section.data(10).logicalSrcIdx = 14;
	  section.data(10).dtTransOffset = 14;
	
	  ;% Control_linear_DWork.Scope2_PWORK.LoggedData
	  section.data(11).logicalSrcIdx = 15;
	  section.data(11).dtTransOffset = 15;
	
	  ;% Control_linear_DWork.Scope3_PWORK.LoggedData
	  section.data(12).logicalSrcIdx = 16;
	  section.data(12).dtTransOffset = 16;
	
	  ;% Control_linear_DWork.Scope4_PWORK.LoggedData
	  section.data(13).logicalSrcIdx = 17;
	  section.data(13).dtTransOffset = 17;
	
	  ;% Control_linear_DWork.ele1_PWORK.LoggedData
	  section.data(14).logicalSrcIdx = 18;
	  section.data(14).dtTransOffset = 18;
	
	  ;% Control_linear_DWork.ele2_PWORK.LoggedData
	  section.data(15).logicalSrcIdx = 19;
	  section.data(15).dtTransOffset = 19;
	
	  ;% Control_linear_DWork.pitch1_PWORK.LoggedData
	  section.data(16).logicalSrcIdx = 20;
	  section.data(16).dtTransOffset = 20;
	
	  ;% Control_linear_DWork.pitch2_PWORK.LoggedData
	  section.data(17).logicalSrcIdx = 21;
	  section.data(17).dtTransOffset = 21;
	
	  ;% Control_linear_DWork.scope1_PWORK.LoggedData
	  section.data(18).logicalSrcIdx = 22;
	  section.data(18).dtTransOffset = 22;
	
	  ;% Control_linear_DWork.scope11_PWORK.LoggedData
	  section.data(19).logicalSrcIdx = 23;
	  section.data(19).dtTransOffset = 23;
	
	  ;% Control_linear_DWork.scope18_PWORK.LoggedData
	  section.data(20).logicalSrcIdx = 24;
	  section.data(20).dtTransOffset = 24;
	
	  ;% Control_linear_DWork.scope2_PWORK.LoggedData
	  section.data(21).logicalSrcIdx = 25;
	  section.data(21).dtTransOffset = 25;
	
	  ;% Control_linear_DWork.scope6_PWORK.LoggedData
	  section.data(22).logicalSrcIdx = 26;
	  section.data(22).dtTransOffset = 26;
	
	  ;% Control_linear_DWork.travel1_PWORK.LoggedData
	  section.data(23).logicalSrcIdx = 27;
	  section.data(23).dtTransOffset = 27;
	
	  ;% Control_linear_DWork.travel2_PWORK.LoggedData
	  section.data(24).logicalSrcIdx = 28;
	  section.data(24).dtTransOffset = 28;
	
	  ;% Control_linear_DWork.TravelAngle1_PWORK.LoggedData
	  section.data(25).logicalSrcIdx = 29;
	  section.data(25).dtTransOffset = 29;
	
	  ;% Control_linear_DWork.TravelAngle2_PWORK.LoggedData
	  section.data(26).logicalSrcIdx = 30;
	  section.data(26).dtTransOffset = 30;
	
	  ;% Control_linear_DWork.TravelAngle1_PWORK_j.LoggedData
	  section.data(27).logicalSrcIdx = 31;
	  section.data(27).dtTransOffset = 31;
	
	  ;% Control_linear_DWork.Scope1_PWORK_p.LoggedData
	  section.data(28).logicalSrcIdx = 32;
	  section.data(28).dtTransOffset = 32;
	
	  ;% Control_linear_DWork.Scope11_PWORK.LoggedData
	  section.data(29).logicalSrcIdx = 33;
	  section.data(29).dtTransOffset = 33;
	
	  ;% Control_linear_DWork.Scope12_PWORK.LoggedData
	  section.data(30).logicalSrcIdx = 34;
	  section.data(30).dtTransOffset = 34;
	
	  ;% Control_linear_DWork.Scope13_PWORK.LoggedData
	  section.data(31).logicalSrcIdx = 35;
	  section.data(31).dtTransOffset = 35;
	
	  ;% Control_linear_DWork.Scope2_PWORK_k.LoggedData
	  section.data(32).logicalSrcIdx = 36;
	  section.data(32).dtTransOffset = 36;
	
	  ;% Control_linear_DWork.Scope3_PWORK_f.LoggedData
	  section.data(33).logicalSrcIdx = 37;
	  section.data(33).dtTransOffset = 37;
	
	  ;% Control_linear_DWork.Scope6_PWORK.LoggedData
	  section.data(34).logicalSrcIdx = 38;
	  section.data(34).dtTransOffset = 38;
	
	  ;% Control_linear_DWork.Scope1_PWORK_n.LoggedData
	  section.data(35).logicalSrcIdx = 39;
	  section.data(35).dtTransOffset = 39;
	
	  ;% Control_linear_DWork.Scope4_PWORK_d.LoggedData
	  section.data(36).logicalSrcIdx = 40;
	  section.data(36).dtTransOffset = 40;
	
	  ;% Control_linear_DWork.Scope7_PWORK.LoggedData
	  section.data(37).logicalSrcIdx = 41;
	  section.data(37).dtTransOffset = 41;
	
	  ;% Control_linear_DWork.Scope9_PWORK.LoggedData
	  section.data(38).logicalSrcIdx = 42;
	  section.data(38).dtTransOffset = 42;
	
	  ;% Control_linear_DWork.scope4_PWORK.LoggedData
	  section.data(39).logicalSrcIdx = 43;
	  section.data(39).dtTransOffset = 43;
	
	  ;% Control_linear_DWork.Scope10_PWORK.LoggedData
	  section.data(40).logicalSrcIdx = 44;
	  section.data(40).dtTransOffset = 44;
	
	  ;% Control_linear_DWork.Scope5_PWORK.LoggedData
	  section.data(41).logicalSrcIdx = 45;
	  section.data(41).dtTransOffset = 45;
	
	  ;% Control_linear_DWork.Scope8_PWORK.LoggedData
	  section.data(42).logicalSrcIdx = 46;
	  section.data(42).dtTransOffset = 46;
	
	  ;% Control_linear_DWork.scope5_PWORK.LoggedData
	  section.data(43).logicalSrcIdx = 47;
	  section.data(43).dtTransOffset = 47;
	
	  ;% Control_linear_DWork.Scope1_PWORK_k.LoggedData
	  section.data(44).logicalSrcIdx = 48;
	  section.data(44).dtTransOffset = 48;
	
	  ;% Control_linear_DWork.FromWorkspace_PWORK_k.TimePtr
	  section.data(45).logicalSrcIdx = 49;
	  section.data(45).dtTransOffset = 49;
	
	  ;% Control_linear_DWork.Scope2_PWORK_f.LoggedData
	  section.data(46).logicalSrcIdx = 50;
	  section.data(46).dtTransOffset = 50;
	
	  ;% Control_linear_DWork.scope11_PWORK_j.LoggedData
	  section.data(47).logicalSrcIdx = 51;
	  section.data(47).dtTransOffset = 51;
	
	  ;% Control_linear_DWork.Scope1_PWORK_d.LoggedData
	  section.data(48).logicalSrcIdx = 52;
	  section.data(48).dtTransOffset = 52;
	
	  ;% Control_linear_DWork.FromWorkspace_PWORK_g.TimePtr
	  section.data(49).logicalSrcIdx = 53;
	  section.data(49).dtTransOffset = 53;
	
	  ;% Control_linear_DWork.Scope2_PWORK_g.LoggedData
	  section.data(50).logicalSrcIdx = 54;
	  section.data(50).dtTransOffset = 54;
	
	  ;% Control_linear_DWork.scope11_PWORK_a.LoggedData
	  section.data(51).logicalSrcIdx = 55;
	  section.data(51).dtTransOffset = 55;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% Control_linear_DWork.sfEvent
	  section.data(1).logicalSrcIdx = 56;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_DWork.sfEvent_b
	  section.data(2).logicalSrcIdx = 57;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(3) = section;
      clear section
      
      section.nData     = 11;
      section.data(11)  = dumData; %prealloc
      
	  ;% Control_linear_DWork.TransportDelay1_IWORK.Tail
	  section.data(1).logicalSrcIdx = 58;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_DWork.TransportDelay_IWORK.Tail
	  section.data(2).logicalSrcIdx = 59;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Control_linear_DWork.FromWorkspace_IWORK.PrevIndex
	  section.data(3).logicalSrcIdx = 60;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Control_linear_DWork.FromWorkspace_IWORK_o.PrevIndex
	  section.data(4).logicalSrcIdx = 61;
	  section.data(4).dtTransOffset = 3;
	
	  ;% Control_linear_DWork.TransportDelay1_IWORK_e.Tail
	  section.data(5).logicalSrcIdx = 62;
	  section.data(5).dtTransOffset = 4;
	
	  ;% Control_linear_DWork.TransportDelay_IWORK_c.Tail
	  section.data(6).logicalSrcIdx = 63;
	  section.data(6).dtTransOffset = 5;
	
	  ;% Control_linear_DWork.FromWorkspace_IWORK_i.PrevIndex
	  section.data(7).logicalSrcIdx = 64;
	  section.data(7).dtTransOffset = 6;
	
	  ;% Control_linear_DWork.TransportDelay2_IWORK.Tail
	  section.data(8).logicalSrcIdx = 65;
	  section.data(8).dtTransOffset = 7;
	
	  ;% Control_linear_DWork.FromWorkspace_IWORK_n.PrevIndex
	  section.data(9).logicalSrcIdx = 66;
	  section.data(9).dtTransOffset = 8;
	
	  ;% Control_linear_DWork.FromWorkspace_IWORK_d.PrevIndex
	  section.data(10).logicalSrcIdx = 67;
	  section.data(10).dtTransOffset = 9;
	
	  ;% Control_linear_DWork.FromWorkspace_IWORK_nh.PrevIndex
	  section.data(11).logicalSrcIdx = 68;
	  section.data(11).dtTransOffset = 10;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(4) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% Control_linear_DWork.is_active_c2_Control_linear
	  section.data(1).logicalSrcIdx = 69;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_DWork.is_active_c1_Control_linear
	  section.data(2).logicalSrcIdx = 70;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(5) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% Control_linear_DWork.isStable
	  section.data(1).logicalSrcIdx = 71;
	  section.data(1).dtTransOffset = 0;
	
	  ;% Control_linear_DWork.doneDoubleBufferReInit
	  section.data(2).logicalSrcIdx = 72;
	  section.data(2).dtTransOffset = 1;
	
	  ;% Control_linear_DWork.isStable_j
	  section.data(3).logicalSrcIdx = 73;
	  section.data(3).dtTransOffset = 2;
	
	  ;% Control_linear_DWork.doneDoubleBufferReInit_m
	  section.data(4).logicalSrcIdx = 74;
	  section.data(4).dtTransOffset = 3;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(6) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (dwork)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    dworkMap.nTotData = nTotData;
    


  ;%
  ;% Add individual maps to base struct.
  ;%

  targMap.paramMap  = paramMap;    
  targMap.signalMap = sigMap;
  targMap.dworkMap  = dworkMap;
  
  ;%
  ;% Add checksums to base struct.
  ;%


  targMap.checksum0 = 2806821005;
  targMap.checksum1 = 1064676375;
  targMap.checksum2 = 4000708734;
  targMap.checksum3 = 2989219468;

