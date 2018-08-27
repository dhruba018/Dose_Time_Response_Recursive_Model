%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Recursive Hybrid Model for Dose-time Drug Response Prediction in Pharamcological Studies
% 'BuildModel' trains the Recursive Hybrid Model for dose-time drug response prediction., 
% 'ModelPredict' predicts the dose-time response using the trained Recursive Hybrid Model at 
%       a specific time point for specific doses available in the training data.
% (c) 2018 S. R. Dhruba, A. Rahman
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef RecursiveHybridModel
    
    properties
        Conc;    X;      Y;
        time_ = 0;      numTree = 100;      rngSeed = 1;         % Default values
        alphaModel;   gammaModel;                                     % RF models for alpha & gamma
    end
    
    methods
        % Class constructor...
        function ParameterModelObject = RecursiveHybridModel(M, X_train, Y_train, varargin)
            if nargin > 0
                nvarargin = numel(varargin);
                inClasses = {class(M), class(X_train), class(Y_train)};
                if (sum(cell2mat(regexpi(inClasses, 'double'))) == nargin - nvarargin)      % Check input types
                    RecursiveHybridModelObject = RecursiveHybridModel;                      % Construct class object
                    RecursiveHybridModelObject.X = X_train;
                    RecursiveHybridModelObject.Y = Y_train;
                    RecursiveHybridModelObject.Conc = M;
                    OptionalFields = fieldnames(RecursiveHybridModelObject);                % Optional inputs
                    for k = 1 : nvarargin
                        eval(['RecursiveHybridModelObject.', OptionalFields{k}]) = varargin{k};
                    end
                    ParameterModelObject = BuildModel(RecursiveHybridModelObject);   % Trained model
                else
                    error('Inputs must be numeric!!!')
                end
            end
        end
        
        % Recursive Hybrid model training... 
        %   Output provides alpha & gamma models...
        function ParameterModelObject = BuildModel(RecursiveHybridModelObject)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FUNCTION
        %      ParameterModelObject = BuildModel(RecursiveHybridModelObject)
        % Trains the Recursive Hybrid Model for dose-time drug response prediction, 
        % where drug response values are predicted at a specific time point using
        % time series data at previous time points.
        % 
        % INPUT ARGUMENTS:
        %      RecursiveHybridModelObject: RecursiveHybridModel object with
        %               the following properties --
        %             Conc:  D x 1 vector of drug concentrations available
        %                       in the training data time series. Prediction will
        %                       be performed at these dose levels.
        %             X:       Training covariate matrix of size (m*D x P x Tx),
        %                       where, m = #subjects, D = #doses, P = #predictors,
        %                       Tx = #covariate time points.
        %             Y:       Training response matrix of size (m*D x Ty),
        %                       where, m = #subjects, D = #doses, 
        %                       Ty = #response time points.
        %             t_:      Previous time point scalar value. By default
        %                       RecursiveHybridModel performs the One-step
        %                       Ahead Prediction i.e., default value is 0.
        %                       [See Dhruba et al. (2018) for more details]
        %        numTree: No. of trees to be used in building the Random Forest 
        %                       (TreeBagger) models. Default value is 100.
        %        rngSeed: Random number generator seed used to reproduce
        %                       the results of RF. Default value is 1.
        %    alphaModel: Model for predicting the Growth parameter [see
        %                       Dhruba et al. (2018)]. A TreeBagger object.
        %  gammaModel: Model for predicting the Scaling parameter [see
        %                       Dhruba et al. (2018)]. A TreeBagger object.
        % 
        % OUTPUT ARGUMENTS:
        %       ParameterModelObject: RecursiveHybridModel object with the
        %               same properties. Used as an input to ModerPredict to 
        %               predict the dose-time response. Passes the following 
        %               properties to ModelPredict for prediction purposes --
        %             Conc: Available drug concentrations in training data.
        %             t_:      Previous time point scalar value to be used
        %                       in prediction.
        %    alphaModel: Model for Growth parameter. A TreeBagger object.
        %  gammaModel: Model for Scaling parameter. A TreeBagger object.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            doses = RecursiveHybridModelObject.Conc;            t_ = RecursiveHybridModelObject.time_;
            x_train = RecursiveHybridModelObject.X;               y_train = RecursiveHybridModelObject.Y;
            n_tree = RecursiveHybridModelObject.numTree;     seed = RecursiveHybridModelObject.rngSeed;
            
            [mD_train, ~, Tx] = size(x_train);            D = numel(doses);       m_train = mD_train / D;   % #training sub
            beta_train = [repmat(doses, [m_train, 1]), (x_train(:, :, Tx) - x_train(:, :, Tx-1))];     % dx/dt @ dt = t - t_
            alpha_train = zeros(mD_train, 1);             gamma_train = zeros(mD_train, 1);
            lb = [0, -2];           ub = [5, 5];                 param0 = [1, 1];       % Optimization parameters
            options = optimoptions(@fmincon, 'display', 'off');                % Optimization options
            for k = 1 : mD_train
                costfun = @(param) RecursiveCostFunction(y_train(k, :), param, t_);
                parameters = fmincon(costfun, param0, [ ], [ ], [ ], [ ], lb, ub, [ ], options);
                alpha_train(k) = parameters(2);            gamma_train(k) = parameters(1);
            end
            rng(seed);      RFalpha = TreeBagger(n_tree, beta_train, alpha_train, 'method', 'regression');
            rng(seed);      RFgamma = TreeBagger(n_tree, beta_train, gamma_train, 'method', 'regression');
            
            % Build output object...
            ParameterModelObject = RecursiveHybridModel;
            ParameterModelObject.Conc = doses;                    ParameterModelObject.time_ = t_;
            ParameterModelObject.X = x_train;                       ParameterModelObject.Y = y_train;
            ParameterModelObject.alphaModel = RFalpha;       ParameterModelObject.gammaModel = RFgamma;
            
            % Cost function definition...
            function CostValue = RecursiveCostFunction(ResponseTimeSeries, RecursivePatameters, t_)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FUNCTION
            %      CostValue = RecursiveCostFunction(ResponseTimeSeries, RecursivePatameters, t_)
            % Cost function definition used in the optimization of time
            % series response data for prediction at further time points.
            % [See Eq. (11) in Dhruba et al. (2018)]
            % 
            % INPUT ARGUMENTS:
            %       ResponseTimeSeries: Ty x 1 vector of Training response
            %                                        time series data.
            %       RecursiveParameters: Recurisve Hybrid Model parameters
            %                                        to be optimized (i.e., alpha & gamma).
            %                               t_:     Previous time point scalar value passed 
            %                                        from BuildModel.
            % 
            % OUTPUT ARGUMENTS:
            %                       CostValue: Cost of optimization. A scalar value.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                y_di = ResponseTimeSeries;                 T = length(ResponseTimeSeries);
                alfa = RecursivePatameters(2);             gama = RecursivePatameters(1);
                J = zeros(T-1, 1);                                 y_di_hat = zeros(T, 1);
                for t = 2 : T               % t = 1 corresponds to time == 0
                    y_di_hat(t) = y_di(t-1) * exp(gama * (1 - exp(-alfa)) * exp(-alfa * t_));       % Estimate from y(t_)
                    J(t-1) = (y_di(t) - y_di_hat(t)).^2;
                end
                CostValue = sum(J);
            end
        end
        
        % Predict response @ t using y(t_)...
        function [RecursiveHybridModelPrediction, AlphaMatrix, GammaMatrix] =...
                                                    ModelPredict(ParameterModelObject, X_test, Y0_test)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FUNCTION
        %      [RecursiveHybridModelPrediction, AlphaMatrix, GammaMatrix] =...
        %                                                     ModelPredict(ParameterModelObject, X_test, Y0_test)
        % Predicts dose-time response using the trained Recursive Hybrid Model 
        % at a specific time point for specific doses available in the training data.
        % 
        % INPUT ARGUMENTS:
        %       ParameterModelObject: RecursiveHybridModel object and the
        %                   output of BuildModel.Passes the following properties for 
        %                   prediction purposes --
        %              Conc: D x 1 vector of available drug concentrations 
        %                       (doses) in training data.
        %                  t_: Previous time point scalar value to be used
        %                       in prediction.
        %    alphaModel: Model for Growth parameter. A TreeBagger object.
        %  gammaModel: Model for Scaling parameter. A TreeBagger object.
        %       X_test: Test covariate matrix of size (m*D x P x Tx),
        %                   where, m = #test subjects, D = #dose, P = #predictors,
        %                   Tx = #covariate time points.
        %      Y0_test: Test response matrix of size (m*D x Ty),
        %                    where, m = #subjects, D = #doses, 
        %                    Ty = #previous response time points.
        %                    Prediction will be performed at t = Ty + 1.
        % 
        % OUTPUT ARGUMENTS:
        %       RecursiveHybridModelPrediction: D x m matrix of Recursive
        %                   Hybrid model prediction. Each column denotes
        %                   the prediction for a single subject at D doses
        %                   and at time = Ty + 1.
        %       AlphaMatrix: D x m matrix of the predicted Growth
        %                   coefficient values.
        %       GammaMatrix: D x m matrix of the predicted Scaling
        %                   coefficient values.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            x_test = X_test;                                                   y0_test = Y0_test;
            doses = ParameterModelObject.Conc;                   t_ = ParameterModelObject.time_;
            RFalpha = ParameterModelObject.alphaModel;      RFgamma = ParameterModelObject.gammaModel;
            
            % Recursion relation...
            Hybrid = @(y_k, alfa, gama, t_k) y_k .* exp(gama .* (1 - exp(-alfa)) .* exp(-alfa * t_k));
            
            % Prediction...
            [mD_test, ~, Tx] = size(x_test);               D = numel(doses);           m_test = mD_test / D;   % #test sub
            beta_test = [repmat(doses, [m_test, 1]), (x_test(:, :, Tx) - x_test(:, :, Tx-1))];     % dx/dt @ dt = t - t_
            alpha_test = predict(RFalpha, beta_test);           alpha_test_hat = alpha_test;
            gamma_test = predict(RFgamma, beta_test);        gamma_test_hat = gamma_test;
            Idx = reshape((1 : mD_test)', [D, m_test]);             % Indices for [D x m_test] format
            y_pred = zeros(mD_test, 1);
            for i = 1 : m_test
                for d = 2 : D
                    % Parameter correction @ d using param(d_)...
                    alpha_test_hat(Idx(d, i)) = max(alpha_test_hat(Idx(d, i)), alpha_test_hat(Idx(d-1, i)));
                    gamma_test_hat(Idx(d, i)) = max(gamma_test_hat(Idx(d, i)), gamma_test_hat(Idx(d-1, i)));
                end
                alfa = alpha_test_hat(Idx(:, i));       gama = gamma_test_hat(Idx(:, i));
                y_pred(Idx(:, i)) = Hybrid(y0_test(Idx(:, i), end), alfa, gama, t_);        % Dose-time prediction
            end
            RecursiveHybridModelPrediction = reshape(y_pred, [D, m_test]);      % One subject per column
            AlphaMatrix = reshape(alpha_test_hat, [D, m_test]);            % Predicted final alpha matrix
            GammaMatrix = reshape(gamma_test_hat, [D, m_test]);         % Predicted final gamma matrix
        end
    end
    
end
