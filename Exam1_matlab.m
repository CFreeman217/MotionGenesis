function [t,VAR,Output] = Exam1_matlab
%===========================================================================
% File: Exam1_matlab.m created May 07 2019 by MotionGenesis 5.9.
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
% % Velocities %%%
% % Partial Velocities %%%
% % Forces %%%
% % Effective Moments %%%
% % Kane's Equations %%%
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
Q1Dt=0; Q2Dt=0; U1Dt=0; U2Dt=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
IAx                             =  0.0;                    % UNITS               Constant
IAy                             =  0.0;                    % UNITS               Constant
IAz                             =  0.0;                    % UNITS               Constant
k                               =  0.0;                    % UNITS               Constant
L                               =  0.0;                    % UNITS               Constant
Ma                              =  0.0;                    % UNITS               Constant
Md                              =  0.0;                    % UNITS               Constant

Q1                              =  0.0;                    % UNITS               Initial Value
Q2                              =  0.0;                    % UNITS               Initial Value
U1                              =  0.0;                    % UNITS               Initial Value
U2                              =  0.0;                    % UNITS               Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  10.0;                   % second              Final Time
tStep                           =  0.05;                   % second              Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
Q1Dt = U1;
Q2Dt = U2;

COEF = zeros( 2, 2 );
COEF(1,1) = (IAx-IAz)*sin(Q2)^2 - 1.49*Ma - IAx - 0.1225*Md*L^2;
COEF(2,2) = -IAy - Ma;
RHS = zeros( 1, 2 );
RHS(1) = 0.7*Ma*U2*(U2+1.428571428571429*cos(Q2)*U1) - 0.142886901662352*k*L^2*sin(Q2)*(-2.449489742783178+L/sqrt(-L^2*(-1.748333333333333+  ...
cos(Q2)))) - 0.8573214099741122*cos(Q2)*(8.009831458900994*Ma-k*L*(-2.449489742783178+L/sqrt(-L^2*(-1.748333333333333+cos(Q2)))))  ...
- 2*(IAx-IAz)*sin(Q2)*cos(Q2)*U1*U2;
RHS(2) = 0.6123724356957945*sin(Q2)*(2*k*L*(-2.449489742783178+L/sqrt(-L^2*(-1.748333333333333+cos(Q2))))-16.01966291780199*Ma-k*L^2*(  ...
-2.449489742783178+L/sqrt(-L^2*(-1.748333333333333+cos(Q2))))) + IAx*sin(Q2)*cos(Q2)*U1^2 - IAz*sin(Q2)*cos(Q2)*U1^2 - 0.7*Ma*U1*(  ...
U2+1.428571428571429*cos(Q2)*U1);
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
      fprintf( 1, '\n Output is in the file Exam1_matlab.1\n' );
      fprintf( 1, '\n Note: To automate plotting, issue the command OutputPlot in MotionGenesis.\n' );
      fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
      fprintf( 1, '    someName = load( ''Exam1_matlab.1'' );\n' );
      fprintf( 1, '    plot( someName(:,1), someName(:,2), ''-'', someName(:,1), someName(:,3), ''--'' )\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             Q1             Q2             U1             U2\n' );
      fprintf( 1,                '%%   (second)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier(1) = fopen('Exam1_matlab.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file Exam1_matlab.1'); end
      fprintf(FileIdentifier(1), '%% FILE: Exam1_matlab.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             Q1             Q2             U1             U2\n' );
      fprintf(FileIdentifier(1), '%%   (second)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
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
end    % End of function Exam1_matlab
%=====================================
