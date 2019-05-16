function [t,VAR,Output] = MGSpinStability3DRigidBodyWithBodyXYZAngles
%===========================================================================
% File: MGSpinStability3DRigidBodyWithBodyXYZAngles.m created May 16 2019 by MotionGenesis 5.9.
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
eventDetectedByIntegratorTerminate1OrContinue0 = [];
q1Dt=0; q2Dt=0; q3Dt=0; wxDt=0; wyDt=0; wzDt=0; theta=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
Ixx                             =  1.381e-6;               % kg*m^2              Constant
Iyy                             =  7.405e-7;               % kg*m^2              Constant
Izz                             =  7.405e-7;               % kg*m^2              Constant

q1                              =  0;                      % degrees             Initial Value
q2                              =  0;                      % degrees             Initial Value
q3                              =  0;                      % degrees             Initial Value
wx                              =  0.2;                    % rad/sec             Initial Value
wy                              =  7.0;                    % rad/sec             Initial Value
wz                              =  0.2;                    % rad/sec             Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  8;                      % sec                 Final Time
tStep                           =  0.01;                   % sec                 Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-7;                 %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;
q1 = q1 * DEGtoRAD;
q2 = q2 * DEGtoRAD;
q3 = q3 * DEGtoRAD;

VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.
PlotOutputFiles;


%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
q1Dt = (wx*cos(q3)-wy*sin(q3))/cos(q2);
q2Dt = wx*sin(q3) + wy*cos(q3);
q3Dt = wz - tan(q2)*(wx*cos(q3)-wy*sin(q3));
wxDt = (Iyy-Izz)*wy*wz/Ixx;
wyDt = -(Ixx-Izz)*wx*wz/Iyy;
wzDt = (Ixx-Iyy)*wx*wy/Izz;

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros( 1, 6 );
VAR(1) = q1;
VAR(2) = q2;
VAR(3) = q3;
VAR(4) = wx;
VAR(5) = wy;
VAR(6) = wz;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
q1 = VAR(1);
q2 = VAR(2);
q3 = VAR(3);
wx = VAR(4);
wy = VAR(5);
wz = VAR(6);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros( 1, 6 );
VARp(1) = q1Dt;
VARp(2) = q2Dt;
VARp(3) = q3Dt;
VARp(4) = wxDt;
VARp(5) = wyDt;
VARp(6) = wzDt;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
theta = acos(cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3));

Output = zeros( 1, 2 );
Output(1) = t;
Output(2) = theta*RADtoDEG;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      fclose( FileIdentifier(1) );
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the file MGSpinStability3DRigidBodyWithBodyXYZAngles.1\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t            theta\n' );
      fprintf( 1,                '%%     (sec)        (degrees)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier(1) = fopen('MGSpinStability3DRigidBodyWithBodyXYZAngles.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file MGSpinStability3DRigidBodyWithBodyXYZAngles.1'); end
      fprintf(FileIdentifier(1), '%% FILE: MGSpinStability3DRigidBodyWithBodyXYZAngles.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t            theta\n' );
      fprintf(FileIdentifier(1), '%%     (sec)        (degrees)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:2) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:2) );  end
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
data = load( 'MGSpinStability3DRigidBodyWithBodyXYZAngles.1' ); 
plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
legend( 'theta (degrees)' );
xlabel('t (sec)');   ylabel('theta (degrees)');   % title('Some plot title');
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


%====================================================================
end    % End of function MGSpinStability3DRigidBodyWithBodyXYZAngles
%====================================================================
