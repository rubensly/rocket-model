% Script to make a basket shot from Rush Rhees library to the Interfaith Chapel (on River Campus, University of Rochester)
close all;clear all;clc

% Set global variables
global lumpedDragCoefficient vehicleMass initialFuelMass initialFuelBurnRate
global expelledMatterVelocity burnTime gravity launchTheta massfuel
 
% Assign variables
vehicleMass = 0.567;        	% kg : mass of basketball ??
initialFuelMass = 0.6;      	% kg
initialFuelBurnRate = .1;   	% kg/s
expelledMatterVelocity = 750;   % m/s
burnTime = 4.07;            	% s
gravity = 9.8;              	% m/s^2
launchTheta = 85;           	% degrees
coefficientOfDrag = 0.47;   	% for sphere
areaOfBluffBody = pi/4*(.749/pi)^2; % m^2
densityOfAir = 1.298;       	% kg/m^3
lumpedDragCoefficient = 1/2*coefficientOfDrag*densityOfAir*areaOfBluffBody;
simTime = 9.5269;           	% s (Simulation time)
 
d_library_chapel = 342.1;   	% meters
radius_ring = .23;          	% meters
 
prompt = {... % these are the prompts shown for the dialog
	'Launch angle (deg)',...
	'Initial fuel mass (kg)',...
	'Initial fuel burn rate (kg/s)',...
	'Time of rocket burn (s)',...
	'Simulation Time (s)',...
	};
dlg_title = 'Input';num_lines = 1; % title and number of lines for inputs
def = {... % these are the default values as strings
	num2str(launchTheta),...	
	num2str(initialFuelMass),...
	num2str(initialFuelBurnRate),...
	num2str(burnTime)...
	num2str(simTime)...
	};
answer = 'x'; % set the starting answer to something notempty
while ~isempty(answer) % while the answer is notempty keep looping
	answer = inputdlg(prompt,dlg_title,num_lines,def); % call the dialog
	if ~isempty(answer) % if the user hits OK, then process it
    	def = answer; % set the current values to the new defaults
    	data = str2double(answer); % convert the strings into numbers
    	% assign to global variables
    	launchTheta = data(1); initialFuelMass = data(2); initialFuelBurnRate = data(3);
    	burnTime = data(4); simTime = data(5) ;
    	
    	[t z] = ode45(@odesys, [0 simTime], [0 0 48.768 0]); % Solve ODE system
    	
    	subplot(2,1,1);plot(z(:,1),z(:,3)); % plot trajectory
    	title('Plot of the basketball trajectory')
    	xlabel('Distance (m)');ylabel('Height (m)');
    	
    	v2 = sqrt(z(:,2).^2+z(:,4).^2); 	
    	subplot(2,1,2);plot(t,v2);      	% plot velocity vs time
    	title('Plot of the basketball velocity')
    	xlabel('Time (s)');ylabel('Velocity (m/s)');
	end
	
	h = z(:,3); % Y-position
	d = z(:,1); % X-position
	% Ball angle at the ring
	theta_end = atand((z(end-1,3)-z(end,3))/(z(end,1)-z(end-1,1)));
	% Terminal velocity
	terminal_vel = sqrt(2*vehicleMass*gravity/densityOfAir/areaOfBluffBody/coefficientOfDrag);
	% Calculate angle as the ball goes through the hoop
	theta_end = atand((z(end-1,3)-z(end,3))/(z(end,1)-z(end-1,1)));
	
	if (d_library_chapel-radius_ring)<=d(end) && d(end)<=(d_library_chapel+radius_ring) && theta_end>=45
    	% Confirm basket if the ball passes through the hoop with an
    	% angle greater than 45 degrees
    	fprintf('\nYes, successful Basketball shot!!!\n')
    	% Display values in table
    	fprintf('%21s %5.4f s\n','Time to basket:',simTime)
    	fprintf('%21s %5.3f kg\n','Fuel burned:',massfuel)
    	fprintf('%21s %5.3f m/s\n','Terminal velocity:',terminal_vel)
    	fprintf('%21s %5.3f m/s\n','Initial velocity:',v2(1))
    	fprintf('%21s %5.3f m\n','Maximum height:',max(z(:,3)))
	else
    	fprintf('Not a basket :( Try again...\n')
	end
end
 
function dz = odesys(t,z)
	global lumpedDragCoefficient vehicleMass initialFuelMass initialFuelBurnRate
	global expelledMatterVelocity burnTime gravity launchTheta massfuel
	
	dz = zeros(4,1);
	if t < burnTime
    	fuelBurnRate = initialFuelBurnRate;
    	massfuel = initialFuelMass - initialFuelBurnRate*t;
    	mass = vehicleMass + massfuel;
	else
    	fuelBurnRate = 0;
    	massfuel = initialFuelMass - initialFuelBurnRate*burnTime;
    	mass = vehicleMass + massfuel;
	end
	if massfuel<0
    	fprintf('Negative mass. Cant calculate\n');
    	%return
	end
	v1 = sqrt(z(2)^2+z(4)^2); % Compute velocity
	if t == 0
    	% CosTheta and SinTheta = value of CosTheta and SinTheta when t=0
    	CosTheta = cos(launchTheta*pi/180);
    	SinTheta = sin(launchTheta*pi/180);
	else
    	v1 = sqrt(z(2)^2+z(4)^2); % Compute velocity
    	% CosTheta and SinTheta = value of CosTheta and SinTheta when t>0
    	CosTheta = z(2)/v1;
    	SinTheta = z(4)/v1;
	end
	% Compute the values for dz to be passed back
	dz(1) = z(2);
	dz(2) = CosTheta*(-lumpedDragCoefficient*v1^2+ expelledMatterVelocity*fuelBurnRate)/mass;
	dz(3) = z(4);
	dz(4) = SinTheta*(-lumpedDragCoefficient.*v1^2+ expelledMatterVelocity*fuelBurnRate)/mass-gravity;
end
