%% Enhanced MATLAB Variable
% -------------------------------------------------------------------------
%  Class v (variable) stores a number with its uncertainty
%  and units. Unit conversion is automatic with other v's.
%  Room for improvement:
%   1. implement matrix compatibility
%   2. automatic unit conversion to simplified units (ie kg*m/s^2 to N)
%   3. upgrade check dimension function to work with multiple dimensions
%  If you add a new non-fundamental unit, make sure that its definition
%  includes all and only the fundamental units that composes it. Also,
%  please do not store large pools of data with v's as they take up much
%  more space and computation time. Adding v(1, 'ft') + v(2, 'ft') takes
%  about 800 times longer than 1 + 2. Note that no steps have been taken to
%  improve the efficiency of v.
%
%  Properties:
%   value - the numerical magnitude of the vector
%   units - a string representing the units of the vector (ft, in, etc)
%   unc - the uncertainty of the value in the same units
%   dimensionless - use this as units for a non-dimensional number
%  Methods:
%   listUnits() - (Static) lists the supported v units
%   dimensionallyEquivalent(v, v) - (Static) returns true/1 or false/0
%   v(number, 'unit string', uncertainty) - constructs a v
%   convertTo('unit string') - returns a new v of equivilant magnitude
%   extract('unit string') - returns a v's value after converting it to
%   the given unit string
%   convertToFundamentals() - returns a new v in fundamental units
%   simplifyUnits() - returns a v with only one unit per dimension
%   checkDimension('dimension character')
%  Examples:
%   v.listUnits() % list supported units
%   x = v(1, 'ft') % make a new v
%   x = x.convertTo('in') %  convert to inches
%   v(1, 'm') + 10 + x % add some items
%   v(1, 'N').convertToFundamentals() % show fundamental units
%   v(1, 'm/ft').simplifyUnits() % simplify units
% 
%  Primary maintainer:
%   Josh Villbrandt (josh@javconcepts.com)
%
%  Revision information:
%   v1.0 (2011-02-06) - completed intial class for USC ADT PlaneTools package
%   v1.1 (2013-02-22) - renamed psv to v and made separate class from PT

%%
classdef v
    %% Variable Properties
    properties
        value;
        units;
        unc = 0;
    end
    
    %% Private Properties
    properties (Access = public)
        unitMap;
    end
    
    %% Constant Properties
    properties (Constant)
        dimensionless = '-';
        % UnitBases is just here for reference.
        unitBases = containers.Map(...
            {'L', 'M', 'T', 'I', 't', 'A', 'i', 'F', 'E', 'p', 'W', 'D',...
            'P', 'r', 'c', 'l', 'q', 'V', 'a'},...
            {'length', 'mass', 'time', 'electric current',...
            'temperature', 'amount of substance', 'luminous intensity',...
            'force*', 'energy*', 'pressure*', 'work*', 'density*',...
            'power*', 'electrical resistance*', 'capacitance*',...
            'inductance*', 'electric charge*', 'electromotive force*',...
            'angular measure^'});
        % Defines a unit's base dimension and its conversion to the SI unit
        % for that base dimension.
        unitDefinitions = containers.Map(...
            {'mm', 'cm', 'm', 'km', 'in', 'ft', 'yd', 'mi', 'nmi', 'ly',...
            'g', 'kg', 'slug', 'lbm', 'oz',...
            'ms', 's', 'min', 'hr', 'd', 'wk', 'mth', 'yr',...
            'mA', 'A',...
            'dC', 'K', 'dF', 'dR',...
            'mol',...
            'cd',...
            'N', 'lb',...
            'J',...
            'Pa',...
            'C', 'mAh', 'Ah',...
            'mV', 'V',...
            'W',...
            'ohm',...
            'rad', 'deg', 'rot',...
            'F'}, {...
            {'L', 0.001},... % mm
            {'L', 0.01},... % cm
            {'L', 1},... % m
            {'L', 1000},... % km
            {'L', 0.0254},... % in
            {'L', 0.3048},... % ft
            {'L', 0.9144},... % yd
            {'L', 1609.344},... % mi
            {'L', 1852},... % nmi
            {'L', 9.461e15},... % ly
            {'M', 0.001},... % g
            {'M', 1},... % kg
            {'M', 14.5939029},... % slug
            {'M', 0.45359237},... % lbm
            {'M', 0.0283495231},... %oz
            {'T', 0.001},... % ms
            {'T', 1},... % s
            {'T', 60},... % min
            {'T', 3600},... % hr
            {'T', 86400},... % d
            {'T', 604800},... % wk
            {'T', 2629743.83},... % mth
            {'T', 31556926},... % yr
            {'I', 31556926},... % mA
            {'I', 31556926},... % A
            {'t', 0},... % dC (requires function)
            {'t', 0},... % K (requires function)
            {'t', 0},... % dF (requires function)
            {'t', 0},... % dR (requires function)
            {'A', 1},... % mol
            {'i', 1},... % cd
            {'F', 1, containers.Map({'kg', 'm', 's'}, {1, 1, -2})},... % N
            {'F', 1, containers.Map({'slug', 'ft', 's'}, {1, 1, -2})},... % lb
            {'E', 1, containers.Map({'kg', 'm', 's'}, {1, 2, -2})},... % J
            {'p', 1, containers.Map({'kg', 'm', 's'}, {1, -1, -2})},... % Pa
            {'q', 1, containers.Map({'A', 's'}, {1, 1})},... % C
            {'q', 3.6, containers.Map({'A', 's'}, {1, 1})},... % mAh
            {'q', 3600, containers.Map({'A', 's'}, {1, 1})},... % Ah
            {'V', 1000, containers.Map({'kg', 'm', 's', 'A'}, {1, 2, -3, -1})},... % mV
            {'V', 1, containers.Map({'kg', 'm', 's', 'A'}, {1, 2, -3, -1})},... % V
            {'P', 1, containers.Map({'kg', 'm', 's'}, {1, 2, -3})},... % W
            {'r', 1, containers.Map({'kg', 'm', 's', 'A'}, {1, 2, -3, -2})},... % Ohm
            {'a', 1},...  % rad
            {'a', (2*pi)/360},... % deg
            {'a', 2*pi},... % rot
            {'c', 1} % F
            });
    end
    
    %% Static Methods
    methods (Static)
        
        %% conversionFactor
        % Calculates conversion factor between complex units. Units should
        % be v's or strings and must be dimensionally equivilant.
        function factor = conversionFactor(vA, vB)
            
            % Convert unit strings to v class
            if(~strcmp(class(vA), class(v)))
                vA = v(1, vA);
            end
            if(~strcmp(class(vB), class(v)))
                vB = v(1, vB);
            end
            
            % Check for special case of identical units
            if(strcmp(vA.units, vB.units))
                factor = 1;
                return;
            end
            
            % Check for special case of one being dimensionless
            ineqStr = ['Units ' vA.units ' and ' vB.units ' are not dimensionally equivalent.'];
            assert(v.dimensionallyEquivalent(vA, vB), 'v:UnequalDims', ineqStr);
            
            % Remove equivalent units
            vAComponents = vA.unitMap;
            vBComponents = vB.unitMap;
            keys = vAComponents.keys();
            for ii = 1 : length(keys)
                if(vBComponents.isKey(keys{ii}))
                    if(vAComponents(keys{ii}) > vBComponents(keys{ii}))
                        vAComponents(keys{ii}) = vAComponents(keys{ii}) - vBComponents(keys{ii});
                        vBComponents.remove(keys{ii});
                    elseif(vAComponents(keys{ii}) < vBComponents(keys{ii}))
                        vBComponents(keys{ii}) = vBComponents(keys{ii}) - vAComponents(keys{ii});
                        vAComponents.remove(keys{ii});
                    else
                        vAComponents.remove(keys{ii});
                        vBComponents.remove(keys{ii});
                    end
                end
            end
            
            % Build conversion factor
            factor = 1;
            keys = vAComponents.keys();
            for ii = 1 : length(keys)
                keysB = vBComponents.keys();
                for jj = 1 : length(keysB)
                    defA = v.unitDefinitions(keys{ii});
                    defB = v.unitDefinitions(keysB{jj});
                    
                    % Check for eq dimensions
                    assert(defA{1} ~= 't', ineqStr); % temperature check
                    if(defA{1} == defB{1})
                        % Conversion factor
                        stepFactor = defA{2} / defB{2};
                        
                        if(vAComponents(keys{ii}) > vBComponents(keysB{jj}))
                            factor = factor * stepFactor^vBComponents(keysB{jj});
                            vAComponents(keys{ii}) = vAComponents(keys{ii}) - vBComponents(keysB{jj});
                            vBComponents.remove(keysB{jj});
                        elseif(vAComponents(keys{ii}) < vBComponents(keysB{jj}))
                            factor = factor * stepFactor^vAComponents(keys{ii});
                            vBComponents(keysB{jj}) = vBComponents(keysB{jj}) - vAComponents(keys{ii});
                            vAComponents.remove(keys{ii});
                        else
                            factor = factor * stepFactor^vAComponents(keys{ii});
                            vAComponents.remove(keys{ii});
                            vBComponents.remove(keysB{jj});
                        end
                    end
                end
            end
            
            % Double check that we caught all units
            assert(isempty(vAComponents) && isempty(vBComponents), 'v:convError', 'Conversion error. Something got screwed up...');
        end % conversionFactor
        
        %% dimensionallyEquivalent
        % return whether two units are dimensionally equivalent or not
        function eq = dimensionallyEquivalent(vA, vB)
            
            % Convert unit strings to v class
            if(~strcmp(class(vA), class(v)))
                % should put a character string assertion here
                vA = v(1, vA);
            end
            if(~strcmp(class(vB), class(v)))
                vB = v(1, vB);
            end
            
            % Check for special case of identical units
            if(strcmp(vA.units, vB.units))
                eq = true;
                return;
            end
            
            % Check for special case of one being dimensionless
            if(strcmp(vA.units, v.dimensionless))
                eq = false;
                return;
            end
            
            % Build dimensions for A, B
            dimensionsA = containers.Map();
            dimensionsB = containers.Map();
            keys = vA.unitMap.keys();
            for ii = 1 : length(keys)
                % Check that unit is known
                assert(v.unitDefinitions.isKey(keys(ii)), ['Unit unknown: ' keys(ii)]);
                
                % Add to dimension
                definition = v.unitDefinitions(keys{ii});
                dimensionString = definition{1};
                if(dimensionsA.isKey(dimensionString))
                    dimensionsA(dimensionString) = dimensionsA(dimensionString) + vA.unitMap(keys{ii});
                else
                    dimensionsA(dimensionString) = vA.unitMap(keys{ii});
                end
            end
            keys = vB.unitMap.keys();
            for ii = 1 : length(keys)
                % Check that unit is known
                assert(v.unitDefinitions.isKey(keys(ii)), ['Unit unknown: ' keys(ii)]);
                
                % Add to dimension
                definition = v.unitDefinitions(keys{ii});
                dimensionString = definition{1};
                if(dimensionsB.isKey(dimensionString))
                    dimensionsB(dimensionString) = dimensionsB(dimensionString) +  + vB.unitMap(keys{ii});
                else
                    dimensionsB(dimensionString) = vB.unitMap(keys{ii});
                end
            end
            
            % Check for equivilent dimmensions
            dimensionsBCopy = dimensionsB;
            dimensions = dimensionsA.keys();
            for ii = 1 : length(dimensions)
                % Check that they both have the same dimension
                if(~dimensionsB.isKey(dimensions{ii}) && dimensionsA(dimensions{ii}) ~= 0)
                    eq = false;
                    return;
                end
                
                % Check that both dimensions are to the same degree
                if(dimensionsA(dimensions{ii}) ~= dimensionsB(dimensions{ii}))
                    eq = false;
                    return;
                end
                
                % Remove from copy
                dimensionsBCopy.remove(dimensions{ii});
            end
            
            % Make sure thereare no remaining dimensions of B
            dimensions = dimensionsBCopy.keys();
            for ii = 1 : length(dimensions)
                if(dimensionsB(dimensions{ii}) > 0)
                    eq = false;
                    return;
                end
            end
                
            eq = true;
            return;
            
        end % dimensionallyEquivalent
        
        %% listUnits
        % Displays currently supported v units
        function listUnits()
            disp('PlaneTools Variable - help()');
            disp('');
            disp('(* denotes derived quantities)');
            disp('([] denotes dimension codes)');
            
            % Loop through fundamental dimensions
            keys = v.unitBases.keys();
            for ii = 1 : length(keys)
                str = [v.unitBases(keys{ii}) ' [' keys{ii} ']:'];
                
                % Loop through units
                defs = v.unitDefinitions.keys();
                for jj = 1 : length(defs)
                    def = v.unitDefinitions(defs{jj});
                    if(def{1} == keys{ii})
                        str = [str ' ' defs{jj}];
                    end
                end
                
                disp(str);
            end
        end % listUnits
    end % methods (Static)
    
    %% Private, static Methods
    methods (Access = private, Static)
        %% strToMap
        % Parses a units sting for the individual units
        function map = strToMap(str)
            map = containers.Map();
            
            if(~isempty(str) && ~strcmp(str, v.dimensionless))
                str = ['*' str];
                startPos = 1;
                while (startPos < length(str)) % at least two characters
                    % Get operation
                    if(str(startPos) == '*')
                        flipOperation = false;
                    elseif(str(startPos) == '/')
                        flipOperation = true;
                    else
                        error('Invalid unit syntax: could not find * or / operator.');
                    end

                    % Exponent
                    if(flipOperation)
                        exponent = -1;
                    else
                        exponent = 1;
                    end
                    
                    % Moving along...
                    startPos = startPos + 1;
                    assert(startPos <= length(str), 'Invalid unit syntax: nothing found after * or / operator.');
                    
                    % Operate recusively if ( is found
                    if(str(startPos) == '(')
                        % Find matching )
                        openParens = strfind(str(startPos:length(str)), '(');
                        closeParens = strfind(str(startPos:length(str)), ')');
                        assert(length(openParens) == length(closeParens), 'Uneven parentheses in unit string.');
                        endParens = 0;
                        for ii = 1 : length(closeParens)
                            if((ii + 1) > length(openParens))
                                endParens = closeParens(length(closeParens));
                            elseif(openParens(ii + 1) > closeParens(ii))
                                endParens = closeParens(ii);
                                break;
                            end
                        end
                        assert(endParens > 0, 'Parentheses algorithm error.');
                        endParens = startPos + endParens - 1;
                        
                        % Generate subStr and check for ^
                        subStr = str((startPos + 1) : (endParens - 1));
                        startPos = endParens + 1;
                        if(startPos < length(str) && str(startPos) == '^')
                            endPos = regexp(str(startPos:length(str)), '[^A-Za-z0-9._()\-\^]');
                            if(isempty(endPos))
                                endPos = length(str);
                            else
                                endPos = startPos - 1 + endPos(1) - 1;
                            end
                            exponent = exponent * str2num(str((startPos+1):endPos));
                            
                            startPos = endPos + 1;
                        end
                        
                        % Process subStr
                        if(~isempty(subStr))
                            % Parse separately and combine
                            subStrMap = v.strToMap(subStr);
                            v.mapTimes(subStrMap, exponent);
                            v.mapPlus(map, subStrMap);
                        end
                    else
                        % Get endPos of unit
                        endPos = regexp(str(startPos:length(str)), '[^A-Za-z0-9._()\-\^]');
                        if(isempty(endPos))
                            endPos = length(str);
                        else
                            endPos = startPos - 1 + endPos(1) - 1;
                        end

                        % Pull out unit
                        carrotPos = startPos - 1 + regexp(str(startPos:endPos), '[\^]');
                        if(~isempty(carrotPos))
                            unit = str(startPos:(carrotPos(1)-1));
                            assert(endPos > carrotPos, 'Invalid unit syntax: nothing found after ^ operator.');
                            exponent = exponent * str2num(str((carrotPos(1)+1):endPos));
                        else
                            unit = str(startPos:endPos);
                        end

                        if(unit(1) ~= '1' && ~strcmp(unit, v.dimensionless))
                            % Check the unit is known
                            assert(v.unitDefinitions.isKey(unit), ['Unit unknown: ' unit]);

                            % Add unit to map
                            if(map.isKey(unit))
                                map(unit) = map(unit) + exponent;
                            else
                                map(unit) = exponent;
                            end
                        end

                        % set new start
                        startPos = endPos + 1;
                    end
                end
            end
        end % strToMap
        
        %% mapToStr
        % Parses a units sting for the individual units
        function str = mapToStr(map)
            posStr = '';
            negStr = '';
            keys = map.keys();
            for ii = 1 : length(keys)
                if(map(keys{ii}) > 0)
                    % Add mult symbol
                    if(~isempty(posStr))
                        posStr = [posStr '*'];
                    end
                    
                    % Add unit
                    posStr = [posStr keys{ii}];
                    
                    % Add power
                    if(map(keys{ii}) ~= 1)
                        posStr = [posStr '^' num2str(map(keys{ii}))];
                    end
                elseif(map(keys{ii}) < 0)
                    % Add unit
                    negStr = [negStr '/' keys{ii}];
                    
                    % Add power
                    if(map(keys{ii}) ~= -1)
                        negStr = [negStr '^' num2str(abs(map(keys{ii})))];
                    end
                end
            end
            
            % Combine strings
            if(isempty(posStr) && isempty(negStr))
                str = v.dimensionless;
            else
                if(isempty(posStr))
                    posStr = '1';
                end
                str = [posStr negStr];
            end
        end % mapToStr
        
        %% mapPlus
        % Add addMap to map
        function mapPlus(map, addMap)
            keys = addMap.keys();
            for ii = 1 : length(keys)
                % add unit powers
                if(map.isKey(keys{ii}))
                    map(keys{ii}) = map(keys{ii}) + addMap(keys{ii});
                else
                    map(keys{ii}) = addMap(keys{ii});
                end
                
                % cleanup
                if(map(keys{ii}) == 0)
                    map.remove(keys{ii});
                end
            end
        end % mapPlus
        
        %% minusMap
        % Subtracts minusMap to map
        function mapMinus(map, minusMap)
            keys = minusMap.keys();
            for ii = 1 : length(keys)
                % add unit powers
                if(map.isKey(keys{ii}))
                    map(keys{ii}) = map(keys{ii}) - minusMap(keys{ii});
                else
                    map(keys{ii}) = -1 * minusMap(keys{ii});
                end
                
                % cleanup
                if(map(keys{ii}) == 0)
                    map.remove(keys{ii});
                end
            end
        end % mapMinus
        
        %% mapTimes
        % multiply map values by b
        function mapTimes(map, b)
            keys = map.keys();
            for ii = 1 : length(keys)
                map(keys{ii}) = map(keys{ii}) * b;
            end
        end % mapTimes
        
        %% duplicateMap
        % Duplicates a maps
        function newMap = duplicateMap(mapA)
            newMap = containers.Map();
            keys = mapA.keys();
            for ii = 1 : length(keys)
                newMap(keys{ii}) = mapA(keys{ii});
            end
        end % duplicateMap
    end % methods (Static)
    
    %% Methods
    methods
        %% Constructor
        function obj = v(value, units, unc)
            if(nargin > 2)
                obj.unc = unc;
            else
                obj.unc = 0;
            end
            
            if(nargin > 1)
                obj.units = units;
            else
                obj.units = v.dimensionless;
            end
            
            if(nargin > 0)
                obj.value = value;
            else
                obj.value = 0;
            end
        end
        
        %% convertTo
        % Convert this v unit to another unit.
        function obj = convertTo(obj, newUnits)
            % Pull out newUnits from v if needed
            if(strcmp(class('newUnits'), 'v'))
                newUnits = newUnits.units;
            end
            % Check for temperature conversion
            if((strcmp(obj.units, 'dC') || strcmp(obj.units, 'K') || strcmp(obj.units, 'dF') || strcmp(obj.units, 'dR')) &&...
                (strcmp(newUnits,'dC') || strcmp(newUnits,'K') || strcmp(newUnits,'dF') || strcmp(newUnits,'dR')))
                % Check to see if the units are different
                if(~strcmp(obj.units, newUnits))
                    if(strcmp(obj.units, 'dC') && strcmp(newUnits, 'K')) obj.value = obj.value + 273.15;
                    elseif(strcmp(obj.units, 'dC') && strcmp(newUnits, 'dF')) obj.value = 9/5*obj.value + 32;
                    elseif(strcmp(obj.units, 'dC') && strcmp(newUnits, 'dR')) obj.value = 9/5*(obj.value + 273.15);
                    elseif(strcmp(obj.units, 'K') && strcmp(newUnits, 'dC')) obj.value = obj.value - 273.15;
                    elseif(strcmp(obj.units, 'K') && strcmp(newUnits, 'dF')) obj.value = 9/5*obj.value - 459.67;
                    elseif(strcmp(obj.units, 'K') && strcmp(newUnits, 'dR')) obj.value = 9/5*obj.value;
                    elseif(strcmp(obj.units, 'dF') && strcmp(newUnits, 'dC')) obj.value = 5/9*(obj.value - 32);
                    elseif(strcmp(obj.units, 'dF') && strcmp(newUnits, 'K')) obj.value = 5/9*(obj.value + 459.67);
                    elseif(strcmp(obj.units, 'dF') && strcmp(newUnits, 'dR')) obj.value = obj.value + 459.67;
                    elseif(strcmp(obj.units, 'dR') && strcmp(newUnits, 'dC')) obj.value = 5/9*(obj.value - 491.67);
                    elseif(strcmp(obj.units, 'dR') && strcmp(newUnits, 'K')) obj.value = 5/9*obj.value;
                    elseif(strcmp(obj.units, 'dR') && strcmp(newUnits, 'dF')) obj.value = obj.value - 459.67;
                    end
                    
                    obj.units = newUnits;
                end
            else
                % Allow for non-fundamental conversions
                objTemp = obj.convertToFundamentals();
                newUnitsTemp = v(1, newUnits);
                newUnitsTemp = newUnitsTemp.convertToFundamentals();
                
                % Simplify units
                objTemp = objTemp.simplifyUnits();
                newUnitsTemp = newUnitsTemp.simplifyUnits();
                
                % Normal conversion
                factor = v.conversionFactor(objTemp.units, newUnitsTemp.units);
                obj.value = objTemp.value / newUnitsTemp.value * factor;
                obj.unc = objTemp.unc / newUnitsTemp.value * factor;
                obj.units = newUnits;
            end
        end % convertTo
        
        %% extract
        % Returns a v's value after converting it to the given unit
        % string
        function value = extract(obj, newUnits)
            obj = obj.convertTo(newUnits);
            value = obj.value;
        end % extract
        
        %% convertToFundamentals
        % Convert this v unit to fundamental unit.
        function obj = convertToFundamentals(obj)
            % Check that this object has units
            if(~isempty(obj.unitMap))
                
                % Change all units to fundamental units if needed
                newUnitMap = v.duplicateMap(obj.unitMap);
                keys = obj.unitMap.keys();
                for ii = 1 : length(keys)
                    assert(v.unitDefinitions.isKey(keys{ii}), ['Unit unknown: ' keys{ii}]);
                    def = v.unitDefinitions(keys{ii});
                    if(length(def) >= 3)
                        % Units
                        addMap = v.duplicateMap(def{3});
                        v.mapTimes(addMap, obj.unitMap(keys{ii}));
                        v.mapPlus(newUnitMap, addMap);
                        newUnitMap.remove(keys{ii});

                        % Value
                        obj.value = def{2} * obj.value;
                        obj.unc = def{2} * obj.unc;
                    end
                end

                % Save new map
                obj.units = newUnitMap;
            end
        end % convertToFundamentals
        
        %% simplifyUnits
        % Combine multiple units of one dimension into one unit.
        function obj = simplifyUnits(obj)
            if(~strcmp(obj.units, v.dimensionless))
                dimensionMap = containers.Map();
                tempMap = v.duplicateMap(obj.unitMap);
                keys = tempMap.keys();
                for ii = 1:length(keys)
                    def = v.unitDefinitions(keys{ii});
                    if(dimensionMap.isKey(def{1}))
                        % This is a different unit in the same dimension, squash it!
                        originalUnit = dimensionMap(def{1});
                        conversionFactor = v.conversionFactor(keys{ii}, originalUnit);
                        conversionFactor = conversionFactor ^ tempMap(keys{ii});
                        obj.value = obj.value * conversionFactor;
                        obj.unc = obj.unc * conversionFactor;
                        
                        % Update UnitMap
                        tempMap(originalUnit) = tempMap(originalUnit) + tempMap(keys{ii});
                        tempMap.remove(keys{ii});
                    else
                        % We haven't seen this unit, lets just save it
                        dimensionMap(def{1}) = keys{ii};
                    end
                end
                
                % Update unit string
                obj.units = v.mapToStr(tempMap);
            end
        end % simplifyUnits
        
        %% Overload - set.value
        % Checks for correct data type.
        function obj = set.value(obj, value)
            assert(~isempty(value) && isnumeric(value), 'v:badValType', 'The value of a v should be numeric.');
            obj.value = value;
        end
        
        %% Overload - set.unc
        % Checks for correct data type.
        function obj = set.unc(obj, value)
            assert(~isempty(value) && isnumeric(value), 'v:badValType', 'The uncertainty of a v should be numeric.');
            obj.unc = value;
        end
        
        %% Overload - set.units
        % Automatic recreating of UnitMap property.
        function obj = set.units(obj, value)
            if(isa(value, 'containers.Map'))
                % Prevents duplicate mapToStr call for internal use
                obj.unitMap = value;
            else
                % Convert string to unit map and save cleaned string
                obj.unitMap = v.strToMap(value);
            end
            
            % Save cleaned string
            obj.units = v.mapToStr(obj.unitMap);
        end
        
        %% Overload - plus
        % Calculates a + b.
        % 
        % The uncertainty propagations is derived as follows:
        % 
        % $$f=a+b$$
        % 
        % $$(\Delta f)^2 = (\frac{\partial f}{\partial a} \Delta a +
        % \frac{\partial f}{\partial b} \Delta b)^2$$
        % 
        % $$(\Delta f)^2 = ( 1 \cdot \Delta a + 1 \cdot \Delta b)^2$$
        %
        % Assuming variables a, b are independant:
        % 
        % $$\Delta f = \sqrt{ (\Delta a)^2 + (\Delta b)^2 }$$
        function obj = plus(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                obj = a;
                b = b.convertTo(a.units);
                obj.value = obj.value + b.value;
                obj.unc = sqrt(a.unc^2 + b.unc^2);
            elseif(strcmp(class(a), class(v)))
                obj = a;
                obj.value = obj.value + b;
            else
                obj = b;
                obj.value = a + obj.value;
            end
        end
        
        %% Overload - minus
        % Calculates a - b.
        % 
        % The uncertainty propagations is derived as follows:
        % 
        % $$f=a-b$$
        % 
        % $$(\Delta f)^2 = (\frac{\partial f}{\partial a} \Delta a +
        % \frac{\partial f}{\partial b} \Delta b)^2$$
        % 
        % $$(\Delta f)^2 = ( 1 \cdot \Delta a - 1 \cdot \Delta b)^2$$
        %
        % Assuming variables a, b are independant:
        % 
        % $$\Delta f = \sqrt{ (\Delta a)^2 + (\Delta b)^2 }$$
        function obj = minus(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                obj = a;
                b = b.convertTo(a.units);
                obj.value = obj.value - b.value;
                obj.unc = sqrt(a.unc^2 + b.unc^2);
            elseif(strcmp(class(a), class(v)))
                obj = a;
                obj.value = obj.value - b;
            else
                obj = b;
                obj.value = a - obj.value;
            end
        end
        
        %% Overload - mtimes
        % Calculates a * b.
        % 
        % The uncertainty propagations is derived as follows:
        % 
        % $$f=a \cdot b$$
        % 
        % $$(\Delta f)^2 = (\frac{\partial f}{\partial a} \Delta a +
        % \frac{\partial f}{\partial b} \Delta b)^2$$
        % 
        % $$(\Delta f)^2 = ( b \Delta a + a \Delta b)^2$$
        %
        % Assuming variables a, b are independant:
        % 
        % $$\Delta f = \sqrt{ (b \Delta a)^2 + (a \Delta b)^2 }$$
        function obj = mtimes(a, b)
        % Overload
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                obj = a;
                obj.units = [a.units '*' b.units];
                obj.value = obj.value * b.value;
                obj.unc = sqrt((b.value*a.unc)^2 + (a.value*b.unc)^2);
            elseif(strcmp(class(a), class(v)))
                obj = a;
                obj.value = obj.value * b;
            else
                obj = b;
                obj.value = a * obj.value;
            end
        end
        
        %% Overload - mrdivide
        % Calculates a / b.
        % 
        % The uncertainty propagations is derived as follows:
        % 
        % $$f=\frac{a}{b}$$
        % 
        % $$(\Delta f)^2 = (\frac{\partial f}{\partial a} \Delta a +
        % \frac{\partial f}{\partial b} \Delta b)^2$$
        % 
        % $$(\Delta f)^2=(\frac{1}{b}\Delta a+\frac{-a}{b^2}\Delta b)^2$$
        %
        % Assuming variables a, b are independant:
        % 
        % $$\Delta f = \sqrt{ (\frac{1}{b}\Delta a)^2 + (\frac{a}{b^2}\Delta b)^2 }$$
        function obj = mrdivide(a, b)
        % Overload
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                obj = a;
                obj.units = [a.units '/(' b.units ')'];
                obj.value = obj.value / b.value;
                obj.unc = sqrt((a.unc/b.value)^2 + (a.value/b.value^2*b.unc)^2);
            elseif(strcmp(class(a), class(v)))
                obj = a;
                obj.value = obj.value / b;
            else
                obj = b;
                obj.units = ['1/(' obj.units ')'];
                obj.value = a / obj.value;
            end
        end
        
        
        %% Overload - mpower
        % Calculates a^b.
        % 
        % The uncertainty propagations is derived as follows:
        % 
        % $$f=a^b$$
        % 
        % $$(\Delta f)^2 = (\frac{\partial f}{\partial a} \Delta a +
        % \frac{\partial f}{\partial b} \Delta b)^2$$
        % 
        % $$(\Delta f)^2=(a^{b-1}b\Delta a + a^b ln(a)\Delta b)^2$$
        %
        % Assuming variables a, b are independant:
        % 
        % $$\Delta f = \sqrt{ ((a^{b-1}b\Delta a)^2 +
        % (a^b ln(a)\Delta b)^2 }$$
        function obj = mpower(a, b)
        % Overload
            if(strcmp(class(b), class(v)))
                assert(b.units == v.dimensionless, 'v:unsupported',...
                    'Taking the power of something with units not yet supported');
            end
            
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                obj = a;
                obj.value = obj.value ^ b.value;
                obj.units = ['(' a.units ')^' num2str(b)];
                obj.unc = sqrt((a.value^(b.value-1)*b.value*a.unc)^2 +...
                    (a.value^b.value*log(a.value)*b.unc)^2);
            elseif(strcmp(class(a), class(v)))
                obj = a;
                obj.value = obj.value ^ b;
                obj.units = ['(' a.units ')^' num2str(b)];
                obj.unc = sqrt((a.value^(b-1)*b*a.unc)^2);
            else
                obj = v();
                obj.value = a ^ b.value; % just normal numbers
                obj.unc = sqrt((a^b.value*log(a)*b.unc)^2);
            end
        end
        
        %% Overload - square root
        % Calculates sqrt(a).
        function obj = sqrt(a)
            obj = mpower(a, 0.5);
        end
        
        %% Overload - less than
        function c = lt(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                b = b.convertTo(a.units);
                c = (a.value + a.unc) < (b.value - b.unc);
            elseif(strcmp(class(a), class(v)))
                c = (a.value + a.unc) < b;
            else
                c = a < (b.value - b.unc);
            end
        end
        
        %% Overload - less than or equal to
        function c = le(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                b = b.convertTo(a.units);
                c = (a.value + a.unc) <= (b.value - b.unc);
            elseif(strcmp(class(a), class(v)))
                c = (a.value + a.unc) <= b;
            else
                c = a <= (b.value - b.unc);
            end
        end
        
        %% Overload - greater than
        function c = gt(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                b = b.convertTo(a.units);
                c = (a.value - a.unc) > (b.value + b.unc);
            elseif(strcmp(class(a), class(v)))
                c = (a.value - a.unc) > b;
            else
                c = a > (b.value + b.unc);
            end
        end
        
        %% Overload - greater than or equal to
        function c = ge(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                b = b.convertTo(a.units);
                c = (a.value + a.unc) >= (b.value - b.unc);
            elseif(strcmp(class(a), class(v)))
                c = (a.value + a.unc) >= b;
            else
                c = a >= (b.value - b.unc);
            end
        end
        
        %% Overload - not equal to
        function c = ne(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                b = b.convertTo(a.units);
                c = ((a.value + a.unc) < (b.value - b.unc)) || ((a.value - a.unc) > (b.value + b.unc));
            elseif(strcmp(class(a), class(v)))
                c = ((a.value + a.unc) < b) || ((a.value - a.unc) > b);
            else
                c = (a < (b.value - b.unc)) || (a > (b.value + b.unc));
            end
        end
        
        %% Overload - equal to
        function c = eq(a, b)
            if(strcmp(class(a), class(v)) && strcmp(class(b), class(v)))
                b = b.convertTo(a.units);
                c = ((a.value + a.unc) >= (b.value - b.unc)) && ((a.value - a.unc) <= (b.value + b.unc));
            elseif(strcmp(class(a), class(v)))
                c = ((a.value + a.unc) >= b) && ((a.value - a.unc) <= b);
            else
                c = (a >= (b.value - b.unc)) && (a <= (b.value + b.unc));
            end
        end
        
        %% Overload - sin
        function Y = sin(X)
            X = X.convertTo('rad');
            Y = v(sin(X.value));
        end
        
        %% Overload - asin
        function Y = asin(X)
            assert(X.checkDimension(v.dimensionless), 'v:BadUnits', 'X must be dimensionless.');
            Y = v(asin(X.value), 'rad');
        end
        
        %% Overload - cos
        function Y = cos(X)
            X = X.convertTo('rad');
            Y = v(cos(X.value));
        end
        
        %% Overload - acos
        function Y = acos(X)
            assert(X.checkDimension(v.dimensionless), 'v:BadUnits', 'X must be dimensionless.');
            Y = v(acos(X.value), 'rad');
        end
        
        %% Overload - tan
        function Y = tan(X)
            X = X.convertTo('rad');
            Y = v(tan(X.value));
        end
        
        %% Overload - atan
        function Y = atan(X)
            assert(X.checkDimension(v.dimensionless), 'v:BadUnits', 'X must be dimensionless.');
            Y = v(atan(X.value), 'rad');
        end
        
        %% Check dimension
        % Check dimension, this is pretty limited right now
        % try x = v(1, 'kg'); x.checkDimension('M')
        function r = checkDimension(obj, dim)
            obj = obj.simplifyUnits();
            
            if(strcmp(obj.units, v.dimensionless))
                r = true;
            else
                assert(obj.unitMap.Count == 1, 'checkDimension currently only works when there is one unit');
                def = v.unitDefinitions(obj.units);
                keys = obj.unitMap.keys;
                if(def{1} == dim && obj.unitMap(keys{1}) == 1) %% second argument recently inserted, not fully tested
                    r = true;
                else
                    r = false;
                end
            end
        end
        
        %% Overload - display
        function str = display(a, digits)
            % Value Display Digits
            if(abs(a.unc) >= eps)
                digits = floor(log10(a.unc));
            elseif(nargin < 2)
                digits = NaN;
            else
%                 digits = floor(log10(a.value));
%                 temp = floor(log10(a.value))
%                 if(temp > digits)
%                     digits = temp;
%                 end
            end
            
            % Value
            if(~isnan(digits))
                if(digits > 0)
                    str = num2str(floor(a.value/digits)*digits);
                else
                    str = sprintf('%6.*f', abs(digits), a.value);
                end
            else
                str = num2str(a.value);
            end
            
            % Units
            if(a.units ~= v.dimensionless)
                str = [str ' ' a.units];
            end
            
            % Unc
            if(abs(a.unc) >= eps)
                if(digits > 0)
                    str = [str ' ± ' num2str(floor(a.unc/digits)*digits)];
                else
                    str = [str ' ± ' sprintf('%6.*f', abs(digits), a.unc)];
                end
            
                % Units
                if(a.units ~= v.dimensionless)
                    str = [str ' ' a.units];
                end
            end
            
            if(nargout < 1)
                disp([inputname(1) ' = ' str]);
            end
        end
    end % methods
end % classdef
