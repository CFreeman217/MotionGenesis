function [t,VAR,Output] = Exam2_matlab
%===========================================================================
% File: Exam2_matlab.m created May 08 2019 by MotionGenesis 5.9.
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
% % Kinematical Equations %%%
% % Rotations %%%
% % Translations %%%
% % Constraints %%%
% % Forces %%%
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
U2=0; Q1Dt=0; Q2Dt=0; Q3Dt=0; U1Dt=0; U3Dt=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
g                               =  9.81;                   % m/s^2               Constant
IA                              =  0.12;                   % kg*m^2              Constant
IC                              =  0.04;                   % kg*m^2              Constant
jB                              =  0.06;                   % kg*m^2              Constant
L                               =  0.8;                    % m                   Constant
Ma                              =  0.5;                    % kg                  Constant
Mb                              =  1.0;                    % kg                  Constant
Mc                              =  0.3;                    % kg                  Constant
rB                              =  0.3;                    % m                   Constant
tB                              =  0.1745;                 % rad                 Constant

Q1                              =  0.2618;                 % rad                 Initial Value
Q2                              =  0.02;                   % m                   Initial Value
Q3                              =  0.2618;                 % rad                 Initial Value
U1                              =  0.0873;                 % rad/s               Initial Value
U3                              =  0.0;                    % rad/s               Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  40.0;                   % second              Final Time
tStep                           =  0.1;                    % second              Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
RADtoDEG = 180.0 / pi;

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
PlotOutputFiles;


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
Q1Dt = U1;
U2 = rB*U3;
Q2Dt = U2;
Q3Dt = U3;

COEF = zeros( 2, 2 );
COEF(1,1) = IC + 0.5625*Mc*L^2;
COEF(1,2) = -0.75*Mc*L*Q2*cos(Q1);
COEF(2,1) = -0.75*Mc*L*Q2*cos(Q1);
COEF(2,2) = IA + jB + Mb*Q2^2 + Ma*(L^2+2*L*rB+Q2^2) + 0.5625*Mc*(1.777777777777778*Q2^2+L^2*(2.666666666666667+sin(Q1))^2)  ...
+ rB*(Ma*rB+Mb*rB+Mc*(rB+1.5*L*(2.666666666666667+sin(Q1))));
RHS = zeros( 1, 2 );
RHS(1) = 0.5625*Mc*L*cos(Q1)*(1.333333333333333*g*cos(tB+Q3)+U3*(2.666666666666667*U2+2.666666666666667*L*U3+L*sin(Q1)*U3));
RHS(2) = rB*U3*(Ma*Q2*U3+Mb*Q2*U3+Mc*(Q2*U3-1.5*L*cos(Q1)*U1)) - 0.75*g*(1.333333333333333*Mb*(rB*sin(tB+Q3)+Q2*cos(tB+Q3))+1.333333333333333*  ...
Ma*(L*sin(tB+Q3)+rB*sin(tB+Q3)+Q2*cos(tB+Q3))+Mc*(1.333333333333333*rB*sin(tB+Q3)+1.333333333333333*Q2*cos(tB+Q3)+L*sin(tB+Q3)*(  ...
2.666666666666667+sin(Q1)))) - 2*Ma*Q2*U2*U3 - 2*Mb*Q2*U2*U3 - 0.75*Mc*(2.666666666666667*Q2*U2*U3+L*Q2*sin(Q1)*U1^2+1.5*L^2*cos(  ...
Q1)*(2.666666666666667+sin(Q1))*U1*U3);
SolutionToAlgebraicEquations = COEF \ transpose(RHS);

% Update variables after uncoupling equations
U1Dt = SolutionToAlgebraicEquations(1);
U3Dt = SolutionToAlgebraicEquations(2);

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 5 );
VAR(1) = Q1;
VAR(2) = Q2;
VAR(3) = Q3;
VAR(4) = U1;
VAR(5) = U3;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
Q1 = VAR(1);
Q2 = VAR(2);
Q3 = VAR(3);
U1 = VAR(4);
U3 = VAR(5);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 5 );
VARp(1) = Q1Dt;
VARp(2) = Q2Dt;
VARp(3) = Q3Dt;
VARp(4) = U1Dt;
VARp(5) = U3Dt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
Output = zeros( 1, 4 );
Output(1) = t;
Output(2) = Q1*RADtoDEG;
Output(3) = Q2;
Output(4) = Q3*RADtoDEG;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      fclose( FileIdentifier(1) );
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the file Exam2_matlab.1\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             Q1             Q2             Q3\n' );
      fprintf( 1,                '%%     (sec)          (deg)           (m)           (deg)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier(1) = fopen('Exam2_matlab.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file Exam2_matlab.1'); end
      fprintf(FileIdentifier(1), '%% FILE: Exam2_matlab.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             Q1             Q2             Q3\n' );
      fprintf(FileIdentifier(1), '%%     (sec)          (deg)           (m)           (deg)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:4) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:4) );  end
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
data = load( 'Exam2_matlab.1' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', 'LineWidth',3 );
legend( 'Q1 (deg)', 'Q2 (m)', 'Q3 (deg)' );
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


%=====================================
end    % End of function Exam2_matlab
%=====================================
