function [t,VAR,Output] = Exam4_matlab
%===========================================================================
% File: Exam4_matlab.m created May 08 2019 by MotionGenesis 5.9.
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
% % Velocities %%%
% % Translations %%%
% % Forces %%%
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
Q1Dt=0; Q2Dt=0; U1Dt=0; U2Dt=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
bA                              =  100;                    % N*s/m               Constant
bB                              =  100;                    % N*s/m               Constant
g                               =  9.81;                   % m/s^2               Constant
k                               =  250;                    % N/m                 Constant
lA                              =  1.5;                    % m                   Constant
lB                              =  2.0;                    % m                   Constant
lko                             =  1.0;                    % m                   Constant
mP                              =  10.0;                   % kg                  Constant
mQ                              =  15.0;                   % kg                  Constant

Q1                              =  1.5708;                 % rad                 Initial Value
Q2                              =  1.5708;                 % rad                 Initial Value
U1                              =  0.0;                    % rad                 Initial Value
U2                              =  0.0;                    % rad                 Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  30;                     % second              Final Time
tStep                           =  0.1;                    % second              Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
PlotOutputFiles;


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
Q1Dt = U1;
Q2Dt = U2;

COEF = zeros( 2, 2 );
COEF(1,1) = mP*lA^2 + mQ*(lA^2+lB^2+2*lA*lB*cos(Q2));
COEF(1,2) = mQ*lB*(lB+lA*cos(Q2));
COEF(2,1) = mQ*lB*(lB+lA*cos(Q2));
COEF(2,2) = mQ*lB^2;
RHS = zeros( 1, 2 );
RHS(1) = -mP*g*lA*sin(Q1) - mQ*g*(lA*sin(Q1)+lB*sin(Q1+Q2)) - bA*U1 - bB*U2 - mQ*lA*lB*sin(Q2)*(U1^2-(U1+U2)^2);
RHS(2) = k*lA*lB*sin(Q2)*(1-lko/sqrt(lA^2+lB^2+2*lA*lB*cos(Q2))) - mQ*g*lB*sin(Q1+Q2) - bB*U2 - mQ*lA*lB*sin(Q2)*U1^2;
SolutionToAlgebraicEquations = COEF \ transpose(RHS);

% Update variables after uncoupling equations
U1Dt = SolutionToAlgebraicEquations(1);
U2Dt = SolutionToAlgebraicEquations(2);

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 4 );
VAR(1) = Q1;
VAR(2) = Q2;
VAR(3) = U1;
VAR(4) = U2;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
Q1 = VAR(1);
Q2 = VAR(2);
U1 = VAR(3);
U2 = VAR(4);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 4 );
VARp(1) = Q1Dt;
VARp(2) = Q2Dt;
VARp(3) = U1Dt;
VARp(4) = U2Dt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
Output = zeros( 1, 5 );
Output(1) = t;
Output(2) = Q1;
Output(3) = Q2;
Output(4) = U1;
Output(5) = U2;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      fclose( FileIdentifier(1) );
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the file Exam4_matlab.1\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             Q1             Q2             U1             U2\n' );
      fprintf( 1,                '%%     (sec)          (rad)          (rad)         (rad/s)        (rad/s)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier(1) = fopen('Exam4_matlab.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file Exam4_matlab.1'); end
      fprintf(FileIdentifier(1), '%% FILE: Exam4_matlab.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             Q1             Q2             U1             U2\n' );
      fprintf(FileIdentifier(1), '%%     (sec)          (rad)          (rad)         (rad/s)        (rad/s)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:5) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:5) );  end
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
data = load( 'Exam4_matlab.1' ); 
plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', data(:,1),data(:,4),'--r', data(:,1),data(:,5),'-m', 'LineWidth',3 );
legend( 'Q1 (rad)', 'Q2 (rad)', 'U1 (rad/s)', 'U2 (rad/s)' );
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
end    % End of function Exam4_matlab
%=====================================
