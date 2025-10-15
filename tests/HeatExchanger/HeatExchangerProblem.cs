namespace MGroup.FEM.Thermal.Tests.HeatExchanger
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.Constitutive.ConvectionDiffusion;
	using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
	using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
	using MGroup.DrugDeliveryModel.Tests.Commons;
	using MGroup.DrugDeliveryModel.Tests.Integration;
	using MGroup.FEM.ConvectionDiffusion.Line;
	using MGroup.MSolve.AnalysisWorkflow;
	using MGroup.MSolve.AnalysisWorkflow.Logging;
	using MGroup.MSolve.AnalysisWorkflow.Transient;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Solution;
	using MGroup.MSolve.Solution.AlgebraicModel;
	using MGroup.MSolve.Solution.LinearSystem;
	using MGroup.NumericalAnalyzers;
	using MGroup.NumericalAnalyzers.Dynamic;
	using MGroup.NumericalAnalyzers.Logging;
	using MGroup.NumericalAnalyzers.Staggered;

	using Xunit;

	public class HeatExchangerProblem
	{
        [Theory]
        [InlineData(40000, 40000, 109, 300, 360, 330, 500)]
        [InlineData(40000, 40000, 109, 360, 15, 100, 3600)]
        [InlineData(40000, 4000, 109, 360, 15, 100, 3600)]
        [InlineData(40000, 15000, 69, 360, 15, 100, 3600)]
        public void DoublePipeHeatExchanger(double heatTransferCoeffTube, double heatTransferCoeffShell, double thermalConductivityWall,
            double egTemp, double waterTemp, double wallTemp, double totalTime)
        {
            #region Data
            // DΑΤΑ
            double heatCapacityCoeffTube = 4085;
			double heatCapacityCoeffShell = 4186;
			double heatCapacityCoeffWall = 380;
			double densityTube = 800;
			double densityShell = 1200;
			double densityWall = 8000;
			double flowRateTube = 1;
			double flowRateShell = -1;
			double length = 1;
			double diameterInnerTube = 0.05;
			double diameterOuterTube = 0.06;
			double diameterInnerShell = 0.1;

			//Initial and boundary conditions
			double thetaShellInitial = egTemp; // shell temperature T(x,0)
			double thetaShellBoundaryL = egTemp; // shell temperature T(L,t)
			double thetaTubeInitial = waterTemp; // tube temperature T(x,0)
			double thetaTubeBoundary0 = waterTemp; // tube temperature T(0,t)
			double thetaWallInitial = wallTemp; // wall temperature T(x,0)
			#endregion

			#region Calculate parameters
			// Heat exchanger PDE parameters
			double k1 = 4 * heatTransferCoeffTube / (densityTube * heatCapacityCoeffTube * diameterInnerTube);
			double k2 = 4 * diameterOuterTube * heatTransferCoeffShell / (densityShell * heatCapacityCoeffShell * (Math.Pow(diameterInnerShell, 2) - Math.Pow(diameterOuterTube, 2)));
			double k3 = 4 * diameterInnerTube * heatTransferCoeffTube / (densityWall * heatCapacityCoeffWall * (Math.Pow(diameterOuterTube, 2) - Math.Pow(diameterInnerTube, 2)));
			double k4 = 4 * diameterOuterTube * heatTransferCoeffShell / (densityWall * heatCapacityCoeffWall * (Math.Pow(diameterOuterTube, 2) - Math.Pow(diameterInnerTube, 2)));

			// Convection - Diffusion - Reaction parameters
			// Equation 1 (Tube)
			double alphaCoeffTube = 1;
			double betaCoeffTube = 0;
			double areaTube = Math.PI * Math.Pow(diameterInnerTube / 2, 2);
			double velocityTube = flowRateTube / (areaTube * densityTube);
			double[] gammaCoeffTube = new double[] { velocityTube };
			double deltaCoeffTube = -k1;
			double[] independentCoeffFactorTube = new double[] { k1 };
			double epsilonCoeffTube = independentCoeffFactorTube[0] * thetaWallInitial;
			IConvectionDiffusionProperties materialTube = new ConvectionDiffusionProperties(alphaCoeffTube, betaCoeffTube, gammaCoeffTube, deltaCoeffTube, epsilonCoeffTube);

			// Equation 2 (Shell)
			double alphaCoeffShell = 1;
			double betaCoeffShell = 0;
			double areaShell = Math.PI * (Math.Pow(diameterInnerShell / 2, 2) - Math.Pow(diameterOuterTube / 2, 2));
			double velocityShell = flowRateShell / (areaShell * densityShell);
			double[] gammaCoeffShell = new double[] { velocityShell };
			double deltaCoeffShell = -k2;
			double[] independentCoeffFactorShell = new double[] { k2 };
			double epsilonCoeffShell = independentCoeffFactorShell[0] * thetaWallInitial;
			IConvectionDiffusionProperties materialShell = new ConvectionDiffusionProperties(alphaCoeffShell, betaCoeffShell, gammaCoeffShell, deltaCoeffShell, epsilonCoeffShell);

			// Equation 3 (Wall)
			double alphaCoeffWall = 1;
			double betaCoeffWall = thermalConductivityWall / (heatCapacityCoeffWall * densityWall);
			double areaWall = Math.PI * (Math.Pow(diameterOuterTube / 2, 2) - Math.Pow(diameterInnerTube / 2, 2));
			double[] gammaCoeffWall = new double[] { 0 };
			double deltaCoeffWall = -(k3 + k4);
			double[] independentCoeffFactorWall = new double[] { k3, k4 };
			double epsilonCoeffWall = (independentCoeffFactorWall[0] * thetaTubeInitial) + (independentCoeffFactorWall[1] * thetaShellInitial);
			IConvectionDiffusionProperties materialWall = new ConvectionDiffusionProperties(alphaCoeffWall, betaCoeffWall, gammaCoeffWall, deltaCoeffWall, epsilonCoeffWall);
			#endregion

			#region Create initial models for tube, shell, wall
			// Create tube model
			var modelTube = new Model();
			modelTube.SubdomainsDictionary.Add(0, new Subdomain(0));

			// Create shell model
			var modelShell = new Model();
			modelShell.SubdomainsDictionary.Add(0, new Subdomain(0));

			// Create wall model
			var modelWall = new Model();
			modelWall.SubdomainsDictionary.Add(0, new Subdomain(0));
			#endregion

			#region Create nodes
			// Create nodes
			int nElements = 10;
			int nNodes = nElements + 1;
			double elementLength = length / (double)nElements;
			List<INode> nodes = new List<INode>();
			int id = 0;
			for (int i = 0; i < nNodes; i++)
			{
				nodes.Add(new Node(id, i * length / nElements));
				modelTube.NodesDictionary.Add(nodes[i].ID, nodes[i]);
				modelShell.NodesDictionary.Add(nodes[i].ID, nodes[i]);
				modelWall.NodesDictionary.Add(nodes[i].ID, nodes[i]);
				id++;
			}
			#endregion

			#region Create tube elements
			// Create tube elements
			id = 0;
			for (int i = 0; i < nElements; i++)
			{
				IConvectionDiffusionElementType element = new ConvectionDiffusionRod(new[] { nodes[i], nodes[i + 1] }, areaTube, materialTube);
				element.ID = id;
				element.SubdomainID = 0;
				modelTube.ElementsDictionary.Add(element.ID, element);
				modelTube.SubdomainsDictionary[element.SubdomainID].Elements.Add(element);
				id++;
			}
			#endregion

			#region Create shell elements
			// Create shell elements
			id = 0;
			for (int i = 0; i < nElements; i++)
			{
				IConvectionDiffusionElementType element = new ConvectionDiffusionRod(new[] { nodes[i], nodes[i + 1] }, areaShell, materialShell);
				element.ID = id;
				element.SubdomainID = 0;
				modelShell.ElementsDictionary.Add(element.ID, element);
				modelShell.SubdomainsDictionary[element.SubdomainID].Elements.Add(element);
				id++;
			}
			#endregion

			#region Create wall elements
			// Create wall elements
			id = 0;
			for (int i = 0; i < nElements; i++)
			{
				IConvectionDiffusionElementType element = new ConvectionDiffusionRod(new[] { nodes[i], nodes[i + 1] }, areaWall, materialWall);
				element.ID = id;
				element.SubdomainID = 0;
				modelWall.ElementsDictionary.Add(element.ID, element);
				modelWall.SubdomainsDictionary[element.SubdomainID].Elements.Add(element);
				id++;
			}
			#endregion

			#region Assign tube conditions
			// Assign tube boundary and initial conditions
			modelTube.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
			   new[]
			   {
					new NodalUnknownVariable(nodes[0], ConvectionDiffusionDof.UnknownVariable, thetaTubeBoundary0),
			   },
			   new INodalConvectionDiffusionNeumannBoundaryCondition[] { }));

			List<NodalInitialUnknownVariable> initialConditionsTube = new List<NodalInitialUnknownVariable>();
			foreach (INode node in modelTube.NodesDictionary.Values)
			{
				if (node != modelTube.NodesDictionary[0])
				{
					initialConditionsTube.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, thetaTubeInitial));
				}

			}
			modelTube.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(initialConditionsTube, new IDomainConvectionDiffusionInitialCondition[] { }));
			#endregion

			#region Assign shell conditions
			// Assign shell boundary and initial conditions
			modelShell.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
			   new[]
			   {
					new NodalUnknownVariable(nodes[nodes.Count - 1], ConvectionDiffusionDof.UnknownVariable, thetaShellBoundaryL),
			   },
			   new INodalConvectionDiffusionNeumannBoundaryCondition[] { }));

			List<NodalInitialUnknownVariable> initialConditionsShell = new List<NodalInitialUnknownVariable>();
			foreach (INode node in modelShell.NodesDictionary.Values)
			{
				if (node != modelShell.NodesDictionary[modelShell.NodesDictionary.Count - 1])
				{
					initialConditionsShell.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, thetaShellInitial));
				}
			}
			modelShell.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(initialConditionsShell, new IDomainConvectionDiffusionInitialCondition[] { }));
			#endregion

			#region Assign wall conditions
			// Assign wall boundary conditions
			//modelWall.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
			//   new[]
			//   {
			//		new NodalUnknownVariable(nodes[0], ConvectionDiffusionDof.UnknownVariable, thetaWallBoundary0),
			//		new NodalUnknownVariable(nodes[nodes.Count - 1], ConvectionDiffusionDof.UnknownVariable, thetaWallBoundaryL),
			//   },
			//   new INodalConvectionDiffusionNeumannBoundaryCondition[] { }));

			List<NodalInitialUnknownVariable> initialConditionsWall = new List<NodalInitialUnknownVariable>();
			foreach (INode node in modelWall.NodesDictionary.Values)
			{
				//if ((node != modelWall.NodesDictionary[0]) && (node != modelWall.NodesDictionary[modelWall.NodesDictionary.Count - 1]))
				//{
				initialConditionsWall.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, thetaWallInitial));
				//}
			}
			modelWall.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(initialConditionsWall, new IDomainConvectionDiffusionInitialCondition[] { }));
			#endregion


			// Coupling

			// Tube
			List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBCTube =
				new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
				{
					(BoundaryAndInitialConditionsUtility.BoundaryConditionCase.LeftDirichlet, new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable }, new double[1][] { new double[1]{ nodes[0].X} }, new double[1] { thetaTubeBoundary0 } ),
				};
			List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBCTube =
				new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();

			// Shell
			List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBCShell =
				new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
				{
					(BoundaryAndInitialConditionsUtility.BoundaryConditionCase.RightDirichlet, new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable }, new double[1][] { new double[1]{ nodes[nodes.Count - 1].X } }, new double[1] { thetaShellBoundaryL } ),
				};
			List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBCShell =
				new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();

			// Wall
			List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBCWall =
				new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
				{
					//(BoundaryAndInitialConditionsUtility.BoundaryConditionCase.LeftDirichlet, new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable }, new double[1][] { new double[3]{ nodes[0].X, nodes[0].Y, nodes[0].Z } }, new double[1] { thetaWallBoundary0 } ),
					//(BoundaryAndInitialConditionsUtility.BoundaryConditionCase.RightDirichlet, new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable }, new double[1][] { new double[3]{ nodes[nodes.Count - 1].X, nodes[nodes.Count - 1].Y, nodes[nodes.Count - 1].Z } }, new double[1] { thetaWallBoundaryL } ),
				};
			List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBCWall =
				new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();

			double timeStep = 10;
			int incrementsPertimeStep = 1;
			int currentTimeStep = 0;

			int nodeIDMonitor = nodes.Count % 2 != 0 ? nodes[((int)nodes.Count / 2) + 1].ID : nodes[(int)nodes.Count / 2].ID;
			double[][] nodeIDValues = new double[3][];
			List<(INode, IDofType)> watchDofs = new List<(INode, IDofType)>();

			nodeIDValues[0] = new double[(int)(totalTime + 1)];
			nodeIDValues[1] = new double[(int)(totalTime + 1)];
			nodeIDValues[2] = new double[(int)(totalTime + 1)];

			var dependedVariableTube = new Dictionary<int, double>();
			var dependedVaraibleShell = new Dictionary<int, double>();
			var dependedVariableWall = new Dictionary<int, double>();
			foreach (var element in modelTube.ElementsDictionary.Values)
			{
				dependedVariableTube.Add(element.ID, thetaTubeInitial);
			}
			foreach (var element in modelShell.ElementsDictionary.Values)
			{
				dependedVaraibleShell.Add(element.ID, thetaShellInitial);
			}
			foreach (var element in modelWall.ElementsDictionary.Values)
			{
				dependedVariableWall.Add(element.ID, thetaWallInitial);
			}

			var independedVariableForTube = new Dictionary<int, double[]>();
			var independedVariableForShell = new Dictionary<int, double[]>();
			var independedVariableForWall = new Dictionary<int, double[]>();
			foreach (var element in modelTube.ElementsDictionary.Values)
			{
				independedVariableForTube.Add(element.ID, new double[] { thetaWallInitial });
			}
			foreach (var element in modelShell.ElementsDictionary.Values)
			{
				independedVariableForShell.Add(element.ID, new double[] { thetaWallInitial });
			}
			foreach (var element in modelWall.ElementsDictionary.Values)
			{
				independedVariableForWall.Add(element.ID, new double[] { thetaTubeInitial, thetaShellInitial });
			}


			var eqModelTube = new EquationModelProviderForStagger(modelTube, nodeIDMonitor, ConvectionDiffusionDof.UnknownVariable,
				alphaCoeffTube, betaCoeffTube, gammaCoeffTube, deltaCoeffTube, independentCoeffFactorTube, convectionDiffusionDirichletBCTube, convectionDiffusionNeumannBCTube, independedVariableForTube, dependedVariableTube);
			var eqModelShell = new EquationModelProviderForStagger(modelShell, nodeIDMonitor, ConvectionDiffusionDof.UnknownVariable,
				alphaCoeffShell, betaCoeffShell, gammaCoeffShell, deltaCoeffShell, independentCoeffFactorShell, convectionDiffusionDirichletBCShell, convectionDiffusionNeumannBCShell, independedVariableForShell, dependedVaraibleShell);
			var eqModelWall = new EquationModelProviderForStagger(modelWall, nodeIDMonitor, ConvectionDiffusionDof.UnknownVariable,
				alphaCoeffWall, betaCoeffWall, gammaCoeffWall, deltaCoeffWall, independentCoeffFactorWall, convectionDiffusionDirichletBCWall, convectionDiffusionNeumannBCWall, independedVariableForWall, dependedVariableWall);

			var equationModel = new CoupledEquationsModel(eqModelTube, eqModelShell, eqModelWall, dependedVariableTube, dependedVaraibleShell,
				dependedVariableWall, timeStep, totalTime, incrementsPertimeStep);
			var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
				equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 20, tolerance: 1E-7);

			List<IAnalysisWorkflowLog> resultsTube = new List<IAnalysisWorkflowLog>();
			List<IAnalysisWorkflowLog> resultsShell = new List<IAnalysisWorkflowLog>();
			List<IAnalysisWorkflowLog> resultsWall = new List<IAnalysisWorkflowLog>();
			for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
			{
				equationModel.CurrentTimeStep = currentTimeStep;
				if (currentTimeStep == 0)
				{
					equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
				}
				staggeredAnalyzer.SolveCurrentStep();

				(equationModel.ParentAnalyzers[0] as NewmarkDynamicAnalyzer).AdvanceStep();
				(equationModel.ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();
				(equationModel.ParentAnalyzers[2] as NewmarkDynamicAnalyzer).AdvanceStep();

				for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
				{
					equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
					equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
				}
				#region logging
				resultsTube.Add(equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]);
				resultsShell.Add(equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]);
				resultsWall.Add(equationModel.ParentAnalyzers[2].ChildAnalyzer.Logs[0]);
				#endregion
			}

			Console.WriteLine();
		}
	}
}
