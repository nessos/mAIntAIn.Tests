using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
using System.Security.AccessControl;
using MGroup.Solvers.AlgebraicModel;
using MGroup.LinearAlgebra.Matrices;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using System.Reflection.PortableExecutable;
using System.Xml.Linq;
using TriangleNet.Tools;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Discretization;
using MGroup.FEM.Thermal.Tests.HeatExchanger;
using MGroup.FEM.ConvectionDiffusion.Line;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class EquationModelProviderForStagger
    {
        /// <summary>
        /// List containing the DIRICHLET boundary conditions for the Convection Diffusion problem.
        /// Item1 : Boundary condition case with respect to the face of the domain (LeftDirichlet, TopDirichlet etc).
        /// Item2 : An StructuralDof array containing the DOFs that are constrained.
        /// Item3 : A jagged array of arrays that contain all the coordinate sets of the bounded dofs. The lenght of
        ///         array is equal to the number of constrained dofs. Each array contains the coordinates of the constrained
        ///         dofs.
        /// Item3 : A double array containing the values of the constrained dofs (1-1 correspondence with the dofs in Item2).
        /// </summary>
        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBC;
        
        /// <summary>
        /// List containing the NEUMANN boundary conditions for the Convection Diffusion problem.
        /// Item1 : Boundary condition case with respect to the face of the domain (RightPointFlux, TopDistributedFlux etc)
        /// Item2 : An StructuralDof array containing the information about the direction of the dofs where the force is
        ///         applied.
        /// Item3 : A jagged array of arrays that contain all the coordinate sets of the bounded dofs. The lenght of
        ///         array is equal to the number of constrained dofs. Each array contains the coordinates of the constrained
        ///         dofs
        /// Item3 : A double array containing the values of the constrained dofs (1-1 correspondence with the dofs in Item2)
        /// </summary>
        private List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBC;

		private Model model;

		private double initialCondition;

		private double capacityCoeff;

		private double diffusionCoeff;

		private double[] convectionCoeff;

		private double[] independentProductionCoeff;

		private double dependentProductionCoeff;

		private double[] independedCoeffFactor;

		public Dictionary<int, double[]> independedVariable { get; set; }

		public Dictionary<int, double> DependedVariable { get; }
		public Model Model { get => model; }

		private int nodeIdToMonitor; //TODO put it where it belongs (coupled7and9eqsSolution.cs)
        
        private ConvectionDiffusionDof dofTypeToMonitor;

        public GlobalAlgebraicModel<Matrix> algebraicModel;

		//TODO: Add thetaP Dictionary
        public EquationModelProviderForStagger(Model model, int nodeIdToMonitor, ConvectionDiffusionDof dofTypeToMonitor,
			double capacityCoeff, double diffusionCoeff, double[] convectionCoeff, double dependentProductionCoeff, double[] independedCoeffFactor,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBC,
            List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBC,
			Dictionary<int, double[]> independedVariable, Dictionary<int, double> dependedVariable)
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-30;

			this.model = model;
			this.convectionDiffusionDirichletBC = convectionDiffusionDirichletBC;
			this.convectionDiffusionNeumannBC = convectionDiffusionNeumannBC;
			this.capacityCoeff = capacityCoeff;
			this.diffusionCoeff = diffusionCoeff;
			this.convectionCoeff = convectionCoeff;
			this.dependentProductionCoeff = dependentProductionCoeff;
			this.independedCoeffFactor = independedCoeffFactor;
			this.independedVariable = independedVariable;
			this.DependedVariable = dependedVariable;
            this.nodeIdToMonitor = nodeIdToMonitor;
            this.dofTypeToMonitor = dofTypeToMonitor;
            
        }

		public Model GetModel()
		{
			//double diffusion = k_th;
			Model updatedModel = new Model();
			updatedModel.SubdomainsDictionary.Add(model.SubdomainsDictionary[0].ID, new Subdomain(0));
			foreach (INode node in model.NodesDictionary.Values)
			{
				updatedModel.NodesDictionary.Add(node.ID, node);
			}

			Dictionary<int, double> independedSource = new Dictionary<int, double>();
			foreach (int key in this.model.ElementsDictionary.Keys)
			{
				double sum = 0;
				for (int i = 0; i < independedCoeffFactor.Length; i++)
				{
					sum += independedCoeffFactor[i] * independedVariable[key][i];
				}
				independedSource.Add(key, sum);
				var elementByType = (ConvectionDiffusionRod)model.ElementsDictionary[key];
				double area = elementByType.CrossSectionArea;
				IConvectionDiffusionProperties material = new ConvectionDiffusionProperties(capacityCoeff, diffusionCoeff, convectionCoeff, dependentProductionCoeff /** DependedVariable[key]*/, independedSource[key]);
				IConvectionDiffusionElementType element = new ConvectionDiffusionRod(new[] { model.ElementsDictionary[key].Nodes[0], model.ElementsDictionary[key].Nodes[1] }, area, material);
				element.ID = model.ElementsDictionary[key].ID;
				element.SubdomainID = model.ElementsDictionary[key].SubdomainID;
				updatedModel.ElementsDictionary.Add(element.ID, element);
				updatedModel.SubdomainsDictionary[element.SubdomainID].Elements.Add(element);
			}
			return updatedModel;
		}

		public void AddBoundaryConditions(Model model)
        {
            BoundaryAndInitialConditionsUtility.AssignConvectionDiffusionDirichletBCsToModel(model, convectionDiffusionDirichletBC, 1e-3);
        }

        public (IParentAnalyzer analyzer, ISolver solver, IChildAnalyzer loadcontrolAnalyzer) GetAppropriateSolverAnalyzerAndLog
        (Model model, double pseudoTimeStep, double pseudoTotalTime, int currentStep, int nIncrements)
        {
            var solverFactory = new DenseMatrixSolver.Factory() { IsMatrixPositiveDefinite = false }; //Dense Matrix Solver solves with zero matrices!
            //var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemConvectionDiffusion(model, algebraicModel);
			var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, provider);
			//var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: 1)
			//{
			//	ResidualTolerance = 1E-8,
			//	MaxIterationsPerIncrement = 100,
			//	NumIterationsForMatrixRebuild = 1
			//};
			//var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

			//loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(new List<(INode node, IDofType dof)>()
			//{(model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor)}, algebraicModel);

			//var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
			var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, linearAnalyzer, timeStep: pseudoTimeStep, totalTime: pseudoTotalTime, false, currentStep: currentStep);
			analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
			var analyzer = analyzerBuilder.Build();
			IList<(INode, IDofType)> watchDofs = new List<(INode, IDofType)>
			{
				(model.NodesDictionary[nodeIdToMonitor], dofTypeToMonitor)
			};
			linearAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);
			//loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs[0], algebraicModel);

			return (analyzer, solver, linearAnalyzer/*loadControlAnalyzer*/);
        }

        public List<double> RetrievePressureSolution(ISolver solver, IChildAnalyzer childAnalyzer, Model model, GlobalAlgebraicModel<Matrix> algebraicModel)
        {
            List<double> p = new List<double>();

            var u = childAnalyzer.CurrentAnalysisResult;

            var currentSolution = algebraicModel.ExtractAllResults(u);

            foreach (var node in model.NodesDictionary.Values)
            {
                p.Add(currentSolution.Data[node.ID, 0]);
            }

            return p;
        }

        internal void UpdatePressureDivergenceDictionary(Dictionary<int, double> updateEpsilon, ISolver solver, IChildAnalyzer childAnalyzer, Model model, GlobalAlgebraicModel<Matrix> algebraicModel)
        {
            var p = childAnalyzer.CurrentAnalysisResult;
			var dirichletBoundaryConditions = algebraicModel.BoundaryConditionsInterpreter.GetDirichletBoundaryConditionsWithNumbering()
				.Select(x => new NodalBoundaryCondition(x.Value.Node, x.Key.DOF, x.Value.Amount));
			foreach (var elem in model.ElementsDictionary.Values)
            {
				var elemSolution = algebraicModel.ExtractElementVector(p, elem);
				elem.MapNodalBoundaryConditionsToElementVector(dirichletBoundaryConditions, elemSolution);
				updateEpsilon[elem.ID] = UpdatePressureAndGradientsOfElement(elemSolution);
			}
           
        }

        private double UpdatePressureAndGradientsOfElement(double[] localDisplacements)
        {
			double solution = 0;
			foreach (double value in localDisplacements)
			{
				solution += value / localDisplacements.Length;
			}

            return solution;
        }

		public void UpdateIndependedVariable(Dictionary<int, double>[] newVariables)
		{
			for (int i = 0; i < newVariables.Length; i++)
			{
				foreach (int key in newVariables[i].Keys)
				{
					this.independedVariable[key][i] = newVariables[i][key];
				}
			}
		}

		public void UpdateDependedVariable(Dictionary<int, double> newVariables)
		{
			foreach (var key in newVariables.Keys)
			{
				DependedVariable[key] = newVariables[key];
			}
		}

		public void SaveStateFromElements(Model model)
		{
			
		}

	}
}
