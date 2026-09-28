      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
      
      PARAMETER(NumElemsWithUvar=51558,maxNumLayers=1)
      PARAMETER(NumWindAngles=51558)
      PARAMETER(maxNumWindAngles = NumElemsWithUvar*maxNumLayers)
      DIMENSION WindAngles(NumElemsWithUvar,maxNumLayers)
      DIMENSION iMaterials(NumElemsWithUvar,maxNumLayers)
      DIMENSION iStressTypes(NumElemsWithUvar,maxNumLayers)
      DIMENSION iFailures(NumElemsWithUvar,maxNumLayers)
      DIMENSION iGlobalElemLabels(NumElemsWithUvar)
      
      COMMON /AngleInfo/WindAngles, iGlobalElemLabels, iMaterials,
     $        iStressTypes, iFailures
      SAVE /AngleInfo/
      
      DIMENSION offsets(1)
      CHARACTER*80 tankNames(1)
      CHARACTER*6 tankName
      CHARACTER*113 fullPathName
      LOGICAL = inorder
C     
      DATA WindAngles/maxNumWindAngles*0.0/
      DATA iMaterials/maxNumWindAngles*0/
      DATA iStressTypes/maxNumWindAngles*0/
      DATA iFailures/maxNumWindAngles*0/
      DATA iGlobalElemLabels/NumElemsWithUvar*0/
C     
      DATA tankNames/"TANK-1"/
      DATA offsets/51558/
C
      IF (LOP .eq. 0) THEN
         fullPathName(1:30) = "F:\abaqus_temp\3D-Hashin-4th-c"
         fullPathName(31:60) = "ylinder-PDM\geometry_zhuning_j"
         fullPathName(61:90) = "oint_rebuild\WCM_Zhuning_Joint"
         fullPathName(91:113) = "_Rebuilt_WindAngles.ang"
         open(unit=101, access='SEQUENTIAL', form='FORMATTED',
     $        status='UNKNOWN', file=fullPathName)
        jrcd = 0
        intNUM = 0
        read(101,*) label, layerNum, iTank, iMat, iStressType, iFailure, windAngle
        kOffset = 1
        numOffset = offsets(kOffset)
        lenName = offsets(kOffset)
        tankName = tankNames(kOffset)
        kElem = 1
        labelPrevious = label
        WindAngles(kElem,layerNum) = windAngle
        iMaterials(kElem,layerNum) = iMat
        iStressTypes(kElem,layerNum) = iStressType
        iFailures(kElem,layerNum) = iFailure
        call GetInternal(tankName, label, 1, intNumber, jrcd) 
        iGlobalElemLabels(kElem) = intNumber
        do j=2,NumWindAngles
           read(101,*) label,layerNum,iTank,iMat,iStressType,iFailure,windAngle
           if (label .eq. labelPrevious) then
              WindAngles(kElem,layerNum) = windAngle
              iMaterials(kElem,layerNum) = iMat
              iStressTypes(kElem,layerNum) = iStressType
              iFailures(kElem,layerNum) = iFailure
           else
              kElem = kElem + 1
              if (j .gt. numOffset) then
                 kOffset = kOffset + 1
                 tankName = tankNames(kOffset)
                 numOffset = offsets(kOffset)
              end if
              WindAngles(kElem,layerNum) = windAngle
              iMaterials(kElem,layerNum) = iMat
              iStressTypes(kElem,layerNum) = iStressType
              iFailures(kElem,layerNum) = iFailure
              call GetInternal(tankName, label, 1, intNumber, jrcd) 
              iGlobalElemLabels(kElem) = intNumber
              labelPrevious = label
           end if
        end do
C
C       Perform bubble sort to spead up lookup in UVARM
C
        inorder = .false.
        do while (inorder .eq. .false.)
          inorder = .true.
          do j = 1, NumElemsWithUvar-1
            if (iGlobalElemLabels(j).gt.iGlobalElemLabels(j+1)) then
              iTemp = iGlobalElemLabels(j)
              iGlobalElemLabels(j) = iGlobalElemLabels(j+1)
              iGlobalElemLabels(j+1) = iTemp
              do k = 1, maxNumLayers
                 xTemp =  WindAngles(j,k)
                 WindAngles(j,k) = WindAngles(j+1,k)
                 WindAngles(j+1,k) = xTemp
                 iTemp = iMaterials(j,k)
                 iMaterials(j,k) = iMaterials(j+1,k)
                 iMaterials(j+1,k) = iTemp
                 iTemp = iStressTypes(j,k)
                 iStressTypes(j,k) = iStressTypes(j+1,k)
                 iStressTypes(j+1,k) = iTemp
                 iTemp = iFailures(j,k)
                 iFailures(j,k) = iFailures(j+1,k)
                 iFailures(j+1,k) = iTemp
              end do
              inorder = .false.
            end if
          end do
        end do
C
        close(101)
C
      END IF
C
      RETURN
      END      
C
C
C
C
      FUNCTION GetInterpolatedValue(NumTemps, PropTemps, NumMaterials, 
     $        maxNumTempsPerMaterial, iMaterial, TEMP, Prop)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION NumTemps(NumMaterials)
      DIMENSION PropTemps(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION Prop(NumMaterials,maxNumTempsPerMaterial)
      REAL*8 Prop_pt

      iPoint = -1
      fraction = 0.0
      numTemperatures = NumTemps(iMaterial)
      if (NumTemps(iMaterial) .eq. 1) then
         iPoint = 1
         fraction = 0.0
      else if (TEMP .le. PropTemps(iMaterial,1)) then
         iPoint = 1
         fraction = 0.0
      else if (TEMP .ge. PropTemps(iMaterial,numTemperatures)) then
         iPoint = numTemperatures
         fraction = 0.0
      else
         tempOld = PropTemps(iMaterial,1)
         k = 2
         do while (iPoint .eq. -1)
            tempNew = PropTemps(iMaterial,k)
            if (TEMP .le. tempNew) then
                iPoint = k-1
                fraction = (TEMP-tempOld)/(tempNew-tempOld)
            else
                tempOld = tempNew 
                k = k + 1
            end if
         end do       
      end if 
C     Perform Interpolation
      if (iPoint .eq. 1 .and. fraction .eq. 0.0) then
          Prop_pt = Prop(iMaterial,1)
      else if (iPoint .eq. numTemperatures) then
          Prop_pt = Prop(iMaterial,numTemperatures)
      else
         curValue = Prop(iMaterial,iPoint) 
         nextValue = Prop(iMaterial,iPoint+1)
         Prop_pt = curValue + fraction*(nextValue-curValue)
      end if
      GetInterpolatedValue = Prop_pt
C
      RETURN
      END FUNCTION
C
C
C
C
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     $            NUVARM,NOEL,NPT,NLAYER,NSPT,KSTEP,KINC,
     $            NDI,NSHR,COORD,JMAC,JMATYP, MATLAYO, LACCFLG)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME,CPNAME
      DIMENSION UVAR(*),TIME(2),DIRECT(3,3),T(3,3),COORD(*),
     $     JMAC(*),JMATYP(*) 
C     USER DEFINED DIMENSION STATEMENTS
      REAL*8 windAngle,cosAngle,sinAngle,c2,s2
      LOGICAL found
      CHARACTER*3 FLGRAY(15)
      DIMENSION AR(15),JAR(15)
      DATA angleConv/57.2957795131/
C     
      PARAMETER(NumElemsWithUvar=51558,maxNumLayers=1)
      PARAMETER(NumWindAngles=51558,NumMaterials=1)
      PARAMETER(maxNumTempsPerMaterial=1)
      DIMENSION WindAngles(NumElemsWithUvar,maxNumLayers)
      DIMENSION iMaterials(NumElemsWithUvar,maxNumLayers)
      DIMENSION iStressTypes(NumElemsWithUvar,maxNumLayers)
      DIMENSION iFailures(NumElemsWithUvar,maxNumLayers)
      DIMENSION iGlobalElemLabels(NumElemsWithUvar)
      COMMON /AngleInfo/WindAngles, iGlobalElemLabels, iMaterials,
     $        iStressTypes, iFailures
C
C     The dimensions of the variables AR and JAR
C     must be set equal to or greater than 15
C
C
C     Declare Elastic Lamina Properties:
C
      DIMENSION NumTemps(NumMaterials)
      DIMENSION MatTemps(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION E11(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION E22(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION E33(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION G12(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION G13(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION G23(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION XNU12(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION XNU13(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION XNU23(NumMaterials,maxNumTempsPerMaterial)
      DIMENSION Alphas(NumMaterials)
      DIMENSION TranShearStrengths(NumMaterials)
C
C     Declare Fail Stress Properties:
C
      PARAMETER(maxNumTempsFailStress=1)
      DIMENSION NumStressTemps(NumMaterials)
      DIMENSION StressTemps(NumMaterials,maxNumTempsFailStress)
      DIMENSION Stress_Xt(NumMaterials,maxNumTempsFailStress)
      DIMENSION Stress_Xc(NumMaterials,maxNumTempsFailStress)
      DIMENSION Stress_Yt(NumMaterials,maxNumTempsFailStress)
      DIMENSION Stress_Yc(NumMaterials,maxNumTempsFailStress)
      DIMENSION Stress_S(NumMaterials,maxNumTempsFailStress)
      DIMENSION Stress_CP(NumMaterials,maxNumTempsFailStress)
      DIMENSION Stress_Biax(NumMaterials,maxNumTempsFailStress)
C
C     Declare Fail Strain Properties:
C
      PARAMETER(maxNumTempsFailStrain=1)
      DIMENSION NumStrainTemps(NumMaterials)
      DIMENSION StrainTemps(NumMaterials,maxNumTempsFailStrain)
      DIMENSION Strain_Xt(NumMaterials,maxNumTempsFailStrain)
      DIMENSION Strain_Xc(NumMaterials,maxNumTempsFailStrain)
      DIMENSION Strain_Yt(NumMaterials,maxNumTempsFailStrain)
      DIMENSION Strain_Yc(NumMaterials,maxNumTempsFailStrain)
      DIMENSION Strain_S(NumMaterials,maxNumTempsFailStrain)
C
C     Begin insert from function writeDeclaration
C
C     Additional user-defined declarations written here
C
C
C     Finished insert from function writeDeclaration
C
      COMMON /Debug/iDebug
      DATA Debug/0/
      DATA NumTemps/1/
      DATA Alphas/1.0/
      DATA TranShearStrengths/1.1000E+02/
      MatTemps = reshape((/
     $ 0.0000E+00
     $ /), shape(MatTemps))
C
C      Fill Material Lamina Properties
C
C      Mat 1) t-700 
C
      E11 = reshape((/
     $ 1.4100E+05
     $ /), shape(E11))
      E22 = reshape((/
     $ 1.1400E+04
     $ /), shape(E22))
      E33 = reshape((/
     $ 1.1400E+04
     $ /), shape(E33))
      G12 = reshape((/
     $ 7.1000E+03
     $ /), shape(G12))
      G13 = reshape((/
     $ 7.1000E+03
     $ /), shape(G13))
      G23 = reshape((/
     $ 7.1000E+03
     $ /), shape(G23))
      XNU12 = reshape((/
     $ 2.8000E-01
     $ /), shape(XNU12))
      XNU13 = reshape((/
     $ 2.8000E-01
     $ /), shape(XNU13))
      XNU23 = reshape((/
     $ 2.8000E-01
     $ /), shape(XNU23))
C
C
C      Fill Fail Stress Properties
C
      DATA NumStressTemps/1/
      StressTemps = reshape((/
     $ 0.0000E+00
     $ /), shape(StressTemps))
C
      Stress_Xt = reshape((/
     $ 2.0800E+03
     $ /), shape(Stress_Xt))
      Stress_Xc = reshape((/
     $ 1.2500E+03
     $ /), shape(Stress_Xc))
      Stress_Yt = reshape((/
     $ 6.0000E+01
     $ /), shape(Stress_Yt))
      Stress_Yc = reshape((/
     $ 2.9000E+02
     $ /), shape(Stress_Yc))
      Stress_S = reshape((/
     $ 1.1000E+02
     $ /), shape(Stress_S))
      Stress_CP = reshape((/
     $ 0.0000E+00
     $ /), shape(Stress_CP))
      Stress_Biax = reshape((/
     $ 0.0000E+00
     $ /), shape(Stress_Biax))
C
C      Fill Fail Strain Properties
C
      DATA NumStrainTemps/1/
      StrainTemps = reshape((/
     $ 0.0000E+00
     $ /), shape(StrainTemps))
C
      Strain_Xt = reshape((/
     $ 1.4750E-02
     $ /), shape(Strain_Xt))
      Strain_Xc = reshape((/
     $ 8.8700E-03
     $ /), shape(Strain_Xc))
      Strain_Yt = reshape((/
     $ 5.2600E-03
     $ /), shape(Strain_Yt))
      Strain_Yc = reshape((/
     $ 2.5440E-02
     $ /), shape(Strain_Yc))
      Strain_S = reshape((/
     $ 1.5490E-02
     $ /), shape(Strain_S))
C     The dimensions of the variables AR and JAR
C     must be set equal to or greater than 15
C
      JRCD = 0
      LOCNUM = 0
      WindAngle = 0.0  
      iMaterial = 0  
      iStressType = 0  
      iFailure = 0  
      iLower = 1
      iUpper = NumElemsWithUvar
      found = .false.
      do while (found .eq. .false.)
         if (iLower .eq. iUpper-1) then
             if (iGlobalElemLabels(iUpper) .eq. NOEL) then
                iMid = iUpper
             else if (iGlobalElemLabels(iLower) .eq. NOEL) then
                iMid = iLower
             else
                return
             end if 
         else
            iMid = (iLower + iUpper)/2
         end if
         if (iGlobalElemLabels(iMid) .eq. NOEL) then
              if (NLAYER .eq. 0) then
                 WindAngle = WindAngles(iMid,1)
                 iMaterial = iMaterials(iMid,1)
                 iStressType = iStressTypes(iMid,1)
                 iFailure = iFailures(iMid,1)
              else
                 WindAngle = WindAngles(iMid,NLAYER)
                 iMaterial = iMaterials(iMid,NLAYER)
                 iStressType = iStressTypes(iMid,NLAYER)
                 iFailure = iFailures(iMid,NLAYER)
              end if
              found = .true.
         else if (iGlobalElemLabels(iMid) .lt. NOEL) then
            iLower = iMid
         else
            iUpper = iMid
         end if
      end do
C
      UVAR(1) = WindAngle*angleConv
      C = COS(WindAngle)
      S = SIN(WindAngle)
      SC = S*C
      C2 = C*C
      S2 = S*S
C
C
C
C     Collect temperature in case of temperature dependence.
      CALL GETVRM('TEMP',AR,JAR,FLGRAY,JRCD,
     $            JMAC,JMATYP,MATLAYO, LACCFLG)
      TEMP = AR(1)
C
C     Collect fiber strains.
C     Returned from GETVRM in global order 11,22,33,12,13,23 (21,31,32 for unsymmetric tensors)
      CALL GETVRM('THE',AR,JAR,FLGRAY,JRCD,
     $            JMAC,JMATYP,MATLAYO, LACCFLG)
      THE11 = AR(1)
      THE22 = AR(2)
      THE33 = AR(3)
      THE23 = AR(6)
      THE13 = AR(5)
      THE12 = AR(4)
      CALL GETVRM('LE',AR,JAR,FLGRAY,JRCD,
     $            JMAC,JMATYP,MATLAYO, LACCFLG)
      EE11 = AR(1) - THE11
      EE22 = AR(2) - THE22
      EE33 = AR(3) - THE33
      EE23 = AR(6) - THE23
      EE13 = AR(5) - THE13
      EE12 = AR(4) - THE12
C     1-dir: Meridional, 2-dir: Hoop, 3-dir: Thru-Thickness
      eps11_pos =  EE11*C2 + EE22*S2 + 2.0*EE12*SC
      eps11_neg =  EE11*C2 + EE22*S2 - 2.0*EE12*SC
      eps22_pos =  EE11*S2 + EE22*C2 - 2.0*EE12*SC
      eps22_neg =  EE11*S2 + EE22*C2 + 2.0*EE12*SC
      eps33     =  EE33
      gamma23_pos =  EE13*S + EE23*C
      gamma23_neg = -EE13*S + EE23*C
      gamma13_pos =  EE13*C + EE23*S
      gamma13_neg =  EE13*C - EE23*S
      gamma12_pos =  -EE11*SC + EE22*SC + EE12*(C2-S2)
      gamma12_neg =   EE11*SC - EE22*SC + EE12*(C2-S2)
C
C     Collect stresses (smeared or lamina based on user selection).
      if (iStressType .eq. 1) then
C        LAMINA Stress Output:
         UVAR(2) = eps11_pos
         UVAR(3) = eps11_neg
         UVAR(4) = eps22_pos
         UVAR(5) = eps22_neg
         UVAR(6) = eps33
         UVAR(7) = gamma23_pos
         UVAR(8) = gamma23_neg
         UVAR(9) = gamma13_pos
         UVAR(10) = gamma13_neg
         UVAR(11) = gamma12_pos
         UVAR(12) = gamma12_neg
C        Collect Material Properties at current Temperature.
         E11_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $            maxNumTempsPerMaterial, iMaterial, TEMP, E11) 
         E22_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $            maxNumTempsPerMaterial, iMaterial, TEMP, E22) 
         E33_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $            maxNumTempsPerMaterial, iMaterial, TEMP, E33) 
         G12_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $            maxNumTempsPerMaterial, iMaterial, TEMP, G12) 
         G13_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $            maxNumTempsPerMaterial, iMaterial, TEMP, G13) 
         G23_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $            maxNumTempsPerMaterial, iMaterial, TEMP, G23) 
         XNU12_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $              maxNumTempsPerMaterial, iMaterial, TEMP, XNU12) 
         XNU13_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $              maxNumTempsPerMaterial, iMaterial, TEMP, XNU13) 
         XNU23_pt = GetInterpolatedValue(NumTemps, MatTemps, NumMaterials, 
     $              maxNumTempsPerMaterial, iMaterial, TEMP, XNU23) 
         XNU21_pt = XNU12_pt*E22_pt/E11_pt 
         XNU31_pt = XNU13_pt*E33_pt/E11_pt 
         XNU32_pt = XNU23_pt*E33_pt/E22_pt 
         S11 =  1./E11_pt
         S12 = -XNU21_pt/E22_pt
         S13 = -XNU31_pt/E33_pt
         S21 = -XNU12_pt/E11_pt
         S22 =  1.0/E22_pt
         S23 = -XNU32_pt/E33_pt
         S31 = -XNU13_pt/E11_pt
         S32 = -XNU23_pt/E22_pt
         S33 =  1.0/E33_pt
         det = S11*(S33*S22-S32*S23)-S21*(S33*S12-s32*S13)+
     $         S31*(S23*S12-s22*S13)
         C11 = (S33*S22 - S32*S23)/det
         C12 = -(S33*S12 - S32*S13)/det
         C13 =  (S23*S12 - S22*S13)/det
         C21 = C12
         C22 = (S11*S33 - S31*S13)/det
         C23 = -(S23*S11 - S21*S13)/det
         C31 = C13
         C32 = C23
         C33 = (S11*S22 - S12*S12)/det
C        Compute lamina stresses
         sig11_pos = C11*eps11_pos + C12*eps22_pos + C13*eps33
         sig11_neg = C11*eps11_neg + C12*eps22_neg + C13*eps33
         sig22_pos = C21*eps11_pos + C22*eps22_pos + C23*eps33
         sig22_neg = C21*eps11_neg + C22*eps22_neg + C23*eps33
         sig33_pos = C31*eps11_pos + C32*eps22_pos + C33*eps33
         sig33_neg = C31*eps11_neg + C32*eps22_neg + C33*eps33
         sig12_pos = gamma12_pos*G12_pt
         sig12_neg = gamma12_neg*G12_pt
         sig13_pos = gamma13_pos*G13_pt
         sig13_neg = gamma13_neg*G13_pt
         sig23_pos = gamma23_pos*G23_pt
         sig23_neg = gamma23_neg*G23_pt
         UVAR(13) = sig11_pos
         UVAR(14) = sig11_neg
         UVAR(15) = sig22_pos
         UVAR(16) = sig22_neg
         UVAR(17) = sig33
         UVAR(18) = sig23_pos
         UVAR(19) = sig23_neg
         UVAR(20) = sig13_pos
         UVAR(21) = sig13_neg
         UVAR(22) = sig12_pos
         UVAR(23) = sig12_neg
C        Check if failure criteria are to be included
         if (iFailure .eq. 1) then
            numFailStressTemps = NumStressTemps(iMaterial)
            if (numFailStressTemps .gt. 0) then
C              Collect *FAIL STRESS material data at current temperature
               Sig_Xt = GetInterpolatedValue(NumStressTemps, StressTemps,
     $           NumMaterials,  maxNumTempsFailStress, iMaterial, TEMP, Stress_Xt)
               Sig_Xc = GetInterpolatedValue(NumStressTemps, StressTemps,
     $           NumMaterials,  maxNumTempsFailStress, iMaterial, TEMP, Stress_Xc)
               Sig_Yt = GetInterpolatedValue(NumStressTemps, StressTemps,
     $           NumMaterials,  maxNumTempsFailStress, iMaterial, TEMP, Stress_Yt)
               Sig_Yc = GetInterpolatedValue(NumStressTemps, StressTemps,
     $           NumMaterials,  maxNumTempsFailStress, iMaterial, TEMP, Stress_Yc)
               Sig_S = GetInterpolatedValue(NumStressTemps, StressTemps,
     $           NumMaterials,  maxNumTempsFailStress, iMaterial, TEMP, Stress_S)
               Sig_CP = GetInterpolatedValue(NumStressTemps, StressTemps,
     $           NumMaterials,  maxNumTempsFailStress, iMaterial, TEMP, Stress_CP)
               Sig_Biax = GetInterpolatedValue(NumStressTemps, StressTemps,
     $           NumMaterials,  maxNumTempsFailStress, iMaterial, TEMP, Stress_Biax)
               if (sig11_pos .gt. 0.0) then
                  X_pos = Sig_Xt
               else
                  X_pos = Sig_Xc
               end if
               if (sig22_pos .gt. 0.0) then
                  Y_pos = Sig_Yt
               else
                  Y_pos = Sig_Yc
               end if
               if (sig11_neg .gt. 0.0) then
                  X_neg = Sig_Xt
               else
                  X_neg = Sig_Xc
               end if
               if (sig22_neg .gt. 0.0) then
                  Y_neg = Sig_Yt
               else
                  Y_neg = Sig_Yc
               end if
               S = Sig_S
               CP = Sig_CP
               Biax = Sig_Biax
               alpha = Alphas(iMat)
               tranShear = TranShearStrengths(iMat)
C              Maximum Stress Theory
               UVAR(24) = max(sig11_pos/X_pos,sig22_pos/Y_pos,abs(sig12_pos/S))
               UVAR(25) = max(sig11_neg/X_neg,sig22_neg/Y_neg,abs(sig12_neg/S))
C              Tsai-Hill Theory
               UVAR(26) = (sig11_pos*sig11_pos)/X_pos**2 - 
     $                    (sig11_pos*sig22_pos)/X_pos**2 + 
     $                    (sig22_pos*sig22_pos)/Y_pos**2 + 
     $                    (sig12_pos*sig12_pos)/S**2
               UVAR(27) = (sig11_neg*sig11_neg)/X_neg**2 - 
     $                    (sig11_neg*sig22_neg)/X_neg**2 + 
     $                    (sig22_neg*sig22_neg)/Y_neg**2 + 
     $                    (sig12_neg*sig12_neg)/S**2
C              Tsai-Wu Theory
               F1 = 1.0/Sig_Xt + 1.0/Sig_Xc
               F2 = 1.0/Sig_Yt + 1.0/Sig_Yc
               F11 = -1.0/(Sig_Xt*Sig_Xc)
               F22 = -1.0/(Sig_Yt*Sig_Yc)
               if (Biax .ne. 0.0) then
                  F12 = .5*(1.0-(F1+F2)*Biax+(F11+F22)*Biax**2)/Biax**2
               else
                  F12 = CP*(F11*F22)**0.5
               endif
               F66 = 1.0/S**2
               UVAR(28) = F1*sig11_pos + F2*sig22_pos + F11*sig11_pos**2 + 
     $                    F22*sig22_pos**2 + F66*sig12_pos**2 + 2.0*F12*sig11_pos*sig22_pos
               UVAR(29) =  F1_neg*sig11_neg + F2*sig22_neg + F11*sig11_neg**2 + 
     $                    F22*sig22_neg**2 + F66*sig12_neg**2 + 2.0*F12*sig11_neg*sig22_neg
C              Azzi-Tsai-Hill Theory
               UVAR(30) = sig11_pos**2/X_pos**2 - abs(sig11_pos*sig22_pos)/X_pos**2 + 
     $                    sig22_pos**2/Y_pos**2 + sig12_pos**2/S**2
               UVAR(31) = sig11_neg**2/X_neg**2 - abs(sig11_neg*sig22_neg)/X_neg**2 + 
     $                    sig22_neg**2/Y_neg**2 + sig12_neg**2/S**2
C              Hashin Fiber Failure Theory
               if (sig11_pos .gt. 0.0) then
                  UVAR(32) = sig11_pos**2/X_pos**2 + F66*alpha*sig12_pos**2
               else
                  UVAR(32) = sig11_pos**2/X_pos**2
               end if
               if (sig11_neg .gt. 0.0) then
                  UVAR(33) = sig11_neg**2/X_neg**2 + F66*alpha*sig12_neg**2
               else
                  UVAR(33) = sig11_neg**2/X_neg**2
               end if
C              Hashin Matrix Failure Theory
               if (sig22_pos .gt. 0.0) then
                  UVAR(34) = sig22_pos**2/Y_pos**2 + F66*alpha*sig12_pos**2
               else
                  UVAR(34) = sig22_pos**2/Y_pos**2
               end if
               if (sig22_neg .gt. 0.0) then
                  UVAR(35) = sig22_neg**2/Y_neg**2 + F66*alpha*sig12_neg**2
               else
                  UVAR(35) = sig22_neg**2/Y_neg**2
               end if
            else
               UVAR(24) = 0.0
               UVAR(25) = 0.0
               UVAR(26) = 0.0
               UVAR(27) = 0.0
               UVAR(28) = 0.0
               UVAR(29) = 0.0
               UVAR(30) = 0.0
               UVAR(31) = 0.0
               UVAR(32) = 0.0
               UVAR(33) = 0.0
               UVAR(34) = 0.0
               UVAR(35) = 0.0
            end if
            numFailStrainTemps = NumStrainTemps(iMaterial)
            if (numFailStrainTemps .gt. 0) then
C              Collect *FAIL STRAIN material data at current temperature
               Eps_Xt = GetInterpolatedValue(NumStrainTemps, StrainTemps,
     $           NumMaterials,  maxNumTempsFailStrain, iMaterial, TEMP, Strain_Xt)
               Eps_Xc = GetInterpolatedValue(NumStrainTemps, StrainTemps,
     $           NumMaterials,  maxNumTempsFailStrain, iMaterial, TEMP, Strain_Xc)
               Eps_Yt = GetInterpolatedValue(NumStrainTemps, StrainTemps,
     $           NumMaterials,  maxNumTempsFailStrain, iMaterial, TEMP, Strain_Yt)
               Eps_Yc = GetInterpolatedValue(NumStrainTemps, StrainTemps,
     $           NumMaterials,  maxNumTempsFailStrain, iMaterial, TEMP, Strain_Yc)
               Eps_S = GetInterpolatedValue(NumStrainTemps, StrainTemps,
     $           NumMaterials,  maxNumTempsFailStrain, iMaterial, TEMP, Strain_S)
               if (eps11_pos .gt. 0.0) then
                  X_pos = Eps_Xt
               else
                  X_pos = Eps_Xc
               end if
               if (eps22_pos .gt. 0.0) then
                  Y_pos = Eps_Yt
               else
                  Y_pos = Eps_Yc
               end if
               if (eps11_neg .gt. 0.0) then
                  X_neg = Eps_Xt
               else
                  X_neg = Eps_Xc
               end if
               if (eps22_neg .gt. 0.0) then
                  Y_neg = Eps_Yt
               else
                  Y_neg = Eps_Yc
               end if
               S = Eps_S
C              Maximum Strain Theory
               UVAR(36) = max(eps11_pos/X_pos,eps22_pos/Y_pos,abs(gamma12_pos/S))
               UVAR(37) = max(eps11_neg/X_neg,eps22_neg/Y_neg,abs(gamma12_neg/S))
            else
               UVAR(36) = 0.0
               UVAR(37) = 0.0
            end if
C           Determine Max of all failure critera
            UVAR(38) = max(UVAR(24),UVAR(26),UVAR(28),UVAR(30),UVAR(32),UVAR(34),UVAR(36))
            UVAR(39) = max(UVAR(25),UVAR(27),UVAR(29),UVAR(31),UVAR(33),UVAR(35),UVAR(37))
         end if
      else
C        Smeared Stress Output:
         CALL GETVRM('S',AR,JAR,FLGRAY,JRCD,
     $        JMAC,JMATYP,MATLAYO, LACCFLG)
         UVAR(2) = eps11_pos
         UVAR(3) = eps11_neg
         UVAR(4) = eps22_pos
         UVAR(5) = eps22_neg
         UVAR(6) = gamma12_pos
         UVAR(7) = gamma12_neg
         UVAR(8) = AR(1)*C2 + AR(2)*S2 + 2.0*AR(4)*SC
         UVAR(9) = AR(1)*C2 + AR(2)*S2 - 2.0*AR(4)*SC
         UVAR(10) = AR(1)*S2 + AR(2)*C2 - 2.0*AR(4)*SC
         UVAR(11) = AR(1)*S2 + AR(2)*C2 + 2.0*AR(4)*SC
         UVAR(12) = -AR(1)*SC + AR(2)*SC - AR(4)*(C2-S2)
         UVAR(13) =  AR(1)*SC - AR(2)*SC - AR(4)*(C2-S2)
      end if
C
C     Begin insert from writeExtra
C
C
C     Additional user-defined code written here
C
C
C     Finished insert from writeExtra
C
      RETURN
      END
