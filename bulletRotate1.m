function [t,VAR,Output] = bulletRotate1
%===========================================================================
% File: bulletRotate1.m created May 16 2019 by MotionGenesis 5.9.
% Portions copyright (c) 2009-2019 Motion Genesis LLC.  Rights reserved.
% MotionGenesis Student Licensee: clayFreeman19. (until February 2022).
% Paid-up MotionGenesis Student licensees are granted the right
% to distribute this code for legal student-academic (non-professional) purposes only,
% provided this copyright notice appears in all copies and distributions.
%===========================================================================
% The software is provided "as is", without warranty of any kind, express or    
% implied, including but not limited to the warranties of merchantability or    
% fitness for a particular purpose. In no event shall the authors, contributors,
% or copyright holders be liable for any claim, damages or other liability,     
% whether in an action of contract, tort, or otherwise, arising from, out of, or
% in connection with the software or the use or other dealings in the software. 
%===========================================================================
% % Environment %%%
% % Projectile %%%
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
BaseRadius=0; density=0; IbearingInline=0; IbearingOffAxis=0; IboatTailInline=0; IbtConeBigInline=0; IbtConeBigOffAxis=0; IbtConeSmallInline=0;
IbtConeSmallOffAxis=0; InertiaBoatTailOffAxis=0; InertiaOgiveOffAxis=0; Inertia_Inline=0; Inertia_OffAxis=0; IogiveInline=0; massBaseConeBig=0;
massBaseConeSmall=0; massBearing=0; massOgive=0; volumeBaseConeBig=0; volumeBaseConeSmall=0; volumeBearing=0; volumeBoatTail=0; volumeOgive=0;
volumeTotal=0; Q1Dt=0; QBwxDt=0; QBwyDt=0; QBwzDt=0; U1Dt=0; UBwxDt=0; UBwyDt=0; UBwzDt=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
boatTailAngle                   =  9.0;                    % deg                 Constant
boatTailLength                  =  0.004445;               % m                   Constant
bulletMass                      =  9.52544e-3;             % kg                  Constant
coeffDrag                       =  0.5;                    % noUnits             Constant
g                               =  9.81;                   % m/s^2               Constant
initialAngle                    =  30;                     % deg                 Constant
oal                             =  0.028956;               % m                   Constant
ogiveLength                     =  0.0194818;              % m                   Constant
r                               =  0.00381;                % m                   Constant

Q1                              =  0;                      % m                   Initial Value
QBwx                            =  0;                      % deg                 Initial Value
QBwy                            =  0;                      % deg                 Initial Value
QBwz                            =  0;                      % deg                 Initial Value
U1                              =  860;                    % m/s                 Initial Value
UBwx                            =  177728.15;              % rad/sec             Initial Value
UBwy                            =  0;                      % rad/sec             Initial Value
UBwz                            =  0;                      % rad/sec             Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  60.0;                   % second              Final Time
tStep                           =  0.01;                   % second              Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
boatTailAngle = boatTailAngle * 0.0174532925199433;
initialAngle = initialAngle * 0.0174532925199433;
QBwx = QBwx * 0.0174532925199433;
QBwy = QBwy * 0.0174532925199433;
QBwz = QBwz * 0.0174532925199433;

% Evaluate constants
BaseRadius = r - boatTailLength/tan(boatTailAngle);
volumeBaseConeSmall = 1.047197551196598*BaseRadius^3*tan(boatTailAngle);
volumeBaseConeBig = 1.047197551196598*r^3*tan(boatTailAngle);
volumeOgive = 1.047197551196598*ogiveLength*r^2;
volumeBoatTail = -1.047197551196598*tan(boatTailAngle)*(BaseRadius^3-r^3);
volumeBearing = pi*r^2*(oal-boatTailLength-ogiveLength);
volumeTotal = volumeBearing + volumeBoatTail + volumeOgive;
density = bulletMass/volumeTotal;
massOgive = density*volumeOgive;
massBearing = density*volumeBearing;
massBaseConeSmall = density*volumeBaseConeSmall;
massBaseConeBig = density*volumeBaseConeBig;
IogiveInline = 0.3*massOgive*r^2;
IbtConeSmallInline = 0.3*massBaseConeSmall*BaseRadius^2;
IbtConeBigInline = 0.3*massBaseConeBig*r^2;
IboatTailInline = IbtConeBigInline - IbtConeSmallInline;
IbearingInline = 0.5*massBearing*r^2;
Inertia_Inline = IbearingInline + IboatTailInline + IogiveInline;
IbearingOffAxis = 0.08333333333333333*massBearing*(3*r^2+(oal-boatTailLength-ogiveLength)^2);
IbtConeSmallOffAxis = 0.1*massBaseConeSmall*BaseRadius^2*(1.5+tan(boatTailAngle)^2);
IbtConeBigOffAxis = 0.1*massBaseConeBig*r^2*(1.5+tan(boatTailAngle)^2);
InertiaBoatTailOffAxis = IbtConeBigOffAxis - IbtConeSmallOffAxis;
InertiaOgiveOffAxis = 0.1*massOgive*(ogiveLength^2+1.5*r^2);
Inertia_OffAxis = IbearingOffAxis + InertiaBoatTailOffAxis + InertiaOgiveOffAxis;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
PlotOutputFiles;


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
Q1Dt = U1;
QBwxDt = UBwx;
QBwyDt = UBwy;
QBwzDt = UBwz;
U1Dt = -g*sin(initialAngle) - coeffDrag*U1/bulletMass;
UBwxDt = -coeffDrag*UBwx/Inertia_Inline;
UBwyDt = (1-Inertia_Inline/Inertia_OffAxis)*UBwx*UBwz;
UBwzDt = (-1+Inertia_Inline/Inertia_OffAxis)*UBwx*UBwy;

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 8 );
VAR(1) = Q1;
VAR(2) = QBwx;
VAR(3) = QBwy;
VAR(4) = QBwz;
VAR(5) = U1;
VAR(6) = UBwx;
VAR(7) = UBwy;
VAR(8) = UBwz;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
Q1 = VAR(1);
QBwx = VAR(2);
QBwy = VAR(3);
QBwz = VAR(4);
U1 = VAR(5);
UBwx = VAR(6);
UBwy = VAR(7);
UBwz = VAR(8);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 8 );
VARp(1) = Q1Dt;
VARp(2) = QBwxDt;
VARp(3) = QBwyDt;
VARp(4) = QBwzDt;
VARp(5) = U1Dt;
VARp(6) = UBwxDt;
VARp(7) = UBwyDt;
VARp(8) = UBwzDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
Output = zeros( 1, 6 );
Output(1) = t;
Output(2) = Q1;

Output(3) = t;
Output(4) = UBwx;
Output(5) = UBwy;
Output(6) = UBwz;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 2 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files bulletRotate1.i  (i=1,2)\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             Q1\n' );
      fprintf( 1,                '%%     (sec)           (m)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros( 1, 2 );
      FileIdentifier(1) = fopen('bulletRotate1.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file bulletRotate1.1'); end
      fprintf(FileIdentifier(1), '%% FILE: bulletRotate1.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             Q1\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)\n\n' );
      FileIdentifier(2) = fopen('bulletRotate1.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file bulletRotate1.2'); end
      fprintf(FileIdentifier(2), '%% FILE: bulletRotate1.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%       t            UBwx           UBwy           UBwz\n' );
      fprintf(FileIdentifier(2), '%%     (sec)        (rad/sec)      (rad/sec)      (rad/sec)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:2) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:2) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(3:6) );  end
end


%===========================================================================
function WriteNumericalData( fileIdentifier, Output )
%===========================================================================
numberOfOutputQuantities = length( Output );
if( numberOfOutputQuantities > 0 ),
   for( i = 1 : numberOfOutputQuantities ),
      fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
   end
   fprintf( fileIdentifier, '\n' );
end
end



%===========================================================================
function PlotOutputFiles
%===========================================================================
if( printIntFile == 0 ),  return;  end
figure;
data = load( 'bulletRotate1.1' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'Q1 (m)' );
xlabel('t (sec)');   ylabel('Q1 (m)');   % title('Some plot title');
clear data;

figure;
data = load( 'bulletRotate1.2' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'UBwx (rad/sec)', 'UBwy (rad/sec)', 'UBwz (rad/sec)' );
xlabel('t (sec)');   % ylabel('Some y-axis label');   title('Some plot title');
clear data;
end



%===========================================================================
function [functionsToEvaluateForEvent, eventTerminatesIntegration1Otherwise0ToContinue, eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1] = EventDetection( t, VAR, uSimulink )
%===========================================================================
% Detects when designated functions are zero or cross zero with positive or negative slope.
% Step 1: Uncomment call to mdlDerivatives and mdlOutputs.
% Step 2: Change functionsToEvaluateForEvent,                      e.g., change  []  to  [t - 5.67]  to stop at t = 5.67.
% Step 3: Change eventTerminatesIntegration1Otherwise0ToContinue,  e.g., change  []  to  [1]  to stop integrating.
% Step 4: Change eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1,  e.g., change  []  to  [1].
% Step 5: Possibly modify function EventDetectedByIntegrator (if eventTerminatesIntegration1Otherwise0ToContinue is 0).
%---------------------------------------------------------------------------
% mdlDerivatives( t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
% mdlOutputs(     t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
functionsToEvaluateForEvent = [];
eventTerminatesIntegration1Otherwise0ToContinue = [];
eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1 = [];
eventDetectedByIntegratorTerminate1OrContinue0 = eventTerminatesIntegration1Otherwise0ToContinue;
end


%===========================================================================
function [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents )
%===========================================================================
isIntegrationFinished = eventDetectedByIntegratorTerminate1OrContinue0( nIndexOfEvents );
if( ~isIntegrationFinished ),
   SetNamedQuantitiesFromMatrix( VAR );
%  Put code here to modify how integration continues.
   VAR = SetMatrixFromNamedQuantities;
end
end



%===========================================================================
function [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile )
%===========================================================================
OdeMatlabOptions = odeset( 'RelTol',relError, 'AbsTol',absError, 'MaxStep',tStep, 'Events',@EventDetection );
t = tInitial;                 epsilonT = 0.001*tStep;                   tFinalMinusEpsilonT = tFinal - epsilonT;
printCounterScreen = 0;       integrateForward = tFinal >= tInitial;    tAtEndOfIntegrationStep = t + tStep;
printCounterFile   = 0;       isIntegrationFinished = 0;
mdlDerivatives( t, VAR, 0 );
while 1,
   if( (integrateForward && t >= tFinalMinusEpsilonT) || (~integrateForward && t <= tFinalMinusEpsilonT) ), isIntegrationFinished = 1;  end
   shouldPrintToScreen = printIntScreen && ( isIntegrationFinished || printCounterScreen <= 0.01 );
   shouldPrintToFile   = printIntFile   && ( isIntegrationFinished || printCounterFile   <= 0.01 );
   if( isIntegrationFinished || shouldPrintToScreen || shouldPrintToFile ),
      Output = mdlOutputs( t, VAR, 0 );
      OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile );
      if( isIntegrationFinished ), break;  end
      if( shouldPrintToScreen ), printCounterScreen = printIntScreen;  end
      if( shouldPrintToFile ),   printCounterFile   = printIntFile;    end
   end
   [TimeOdeArray, VarOdeArray, timeEventOccurredInIntegrationStep, nStatesArraysAtEvent, nIndexOfEvents] = ode45( @mdlDerivatives, [t tAtEndOfIntegrationStep], VAR, OdeMatlabOptions, 0 );
   if( isempty(timeEventOccurredInIntegrationStep) ),
      lastIndex = length( TimeOdeArray );
      t = TimeOdeArray( lastIndex );
      VAR = VarOdeArray( lastIndex, : );
      printCounterScreen = printCounterScreen - 1;
      printCounterFile   = printCounterFile   - 1;
      if( abs(tAtEndOfIntegrationStep - t) >= abs(epsilonT) ), warning('numerical integration failed'); break;  end
      tAtEndOfIntegrationStep = t + tStep;
      if( (integrateForward && tAtEndOfIntegrationStep > tFinal) || (~integrateForward && tAtEndOfIntegrationStep < tFinal) ) tAtEndOfIntegrationStep = tFinal;  end
   else
      t = timeEventOccurredInIntegrationStep( 1 );    % time  at firstEvent = 1 during this integration step.
      VAR = nStatesArraysAtEvent( 1, : );             % state at firstEvent = 1 during this integration step.
      printCounterScreen = 0;
      printCounterFile   = 0;
      [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents(1) );
   end
end
end


%======================================
end    % End of function bulletRotate1
%======================================
