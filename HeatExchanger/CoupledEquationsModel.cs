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
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class CoupledEquationsModel
   {
       public EquationModelProviderForStagger Equation1 { get; set; }
       public EquationModelProviderForStagger Equation2 { get; set; }
	   public EquationModelProviderForStagger Equation3 { get; set; }

		public Model[] model;



       //TODO put analysis time stepping where it belongs1 (pithanws sto Coupled7and9eqsSolution.cs h Coupled7and9eqsModel.cs)
       private GenericAnalyzerState[] analyzerStates, nlAnalyzerStates;
       private IParentAnalyzer[] parentAnalyzers;
       private IChildAnalyzer[] nlAnalyzers;
       private ISolver[] parentSolvers;


       public int CurrentTimeStep { get; set; }

       public GenericAnalyzerState[] AnalyzerStates => analyzerStates;
       public GenericAnalyzerState[] NLAnalyzerStates => nlAnalyzerStates;
       public IParentAnalyzer[] ParentAnalyzers => parentAnalyzers;
       public IChildAnalyzer[] NLAnalyzers => nlAnalyzers;
       public ISolver[] ParentSolvers => parentSolvers;

       private Dictionary<int, double> lambda;
       private Dictionary<int, double> thetaTube;
       private Dictionary<int, double> thetaShell;
	   private Dictionary<int, double> thetaWall;

	   private Dictionary<int, double[][]> SolidVelocityAtElementGaussPoints = new Dictionary<int, double[][]>();
       
       private double timeStep;
       private double totalTime;

       private int incrementsPerStep;

       public CoupledEquationsModel(EquationModelProviderForStagger equation1,
	   EquationModelProviderForStagger equation2, EquationModelProviderForStagger equation3,
		   Dictionary<int, double> epsilon1, Dictionary<int, double> epsilon2, Dictionary<int, double> epsilon3,
		   double timeStep, double totalTime, int incrementsPerStep)
       {
			Equation1 = equation1;
			Equation2 = equation2;
			Equation3 = equation3;
			IsoparametricJacobian3D.DeterminantTolerance = 1e-20;


           analyzerStates = new GenericAnalyzerState[3];
           nlAnalyzerStates = new GenericAnalyzerState[3];
           parentAnalyzers = new IParentAnalyzer[3];
           nlAnalyzers = new IChildAnalyzer[3];
           parentSolvers = new ISolver[3];

		   this.thetaTube = epsilon1;
		   this.thetaShell = epsilon2;
		   this.thetaWall = epsilon3;

		   this.timeStep = timeStep;
           this.totalTime = totalTime;
           this.incrementsPerStep = incrementsPerStep;

           // intialize array ofm models1.
           model = new Model[3];
       }


       public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
       {
            Equation1.UpdatePressureDivergenceDictionary(thetaTube, ParentSolvers[0], NLAnalyzers[0], model[0], Equation1.algebraicModel);
			Equation2.UpdatePressureDivergenceDictionary(thetaShell, ParentSolvers[1], NLAnalyzers[1], model[1], Equation2.algebraicModel);
			Equation3.UpdatePressureDivergenceDictionary(thetaWall, ParentSolvers[2], NLAnalyzers[2], model[2], Equation3.algebraicModel);
			
			Equation1.UpdateIndependedVariable(new Dictionary<int, double>[] { thetaWall });
			Equation2.UpdateIndependedVariable(new Dictionary<int, double>[] { thetaWall });
			Equation3.UpdateIndependedVariable(new Dictionary<int, double>[] { thetaTube, thetaShell });

			//Equation1.UpdateDependedVariable(thetaTube);
			//Equation2.UpdateDependedVariable(thetaShell);
			//Equation3.UpdateDependedVariable(thetaWall);

			(ParentAnalyzers[0] as NewmarkDynamicAnalyzer).AdvanceStep();
			(ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();
			(ParentAnalyzers[2] as NewmarkDynamicAnalyzer).AdvanceStep();

			model = new Model[3];

			//Create model for eq78 (fluid pressure)
			model[0] = Equation1.GetModel();
			Equation1.AddBoundaryConditions(model[0]);
			if (CurrentTimeStep == 0)
			{
				model[0].InitialConditions.Add(Equation1.Model.InitialConditions[0]);
			}

            (analyzers[0], solvers[0], nlAnalyzers[0]) = Equation1.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);
			
			model[1] = Equation2.GetModel();
			Equation2.AddBoundaryConditions(model[1]);
			if (CurrentTimeStep == 0)
			{
				model[1].InitialConditions.Add(Equation2.Model.InitialConditions[0]);
			}
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Equation2.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

			model[2] = Equation3.GetModel();
			Equation3.AddBoundaryConditions(model[2]);
			if (CurrentTimeStep == 0)
			{
				model[2].InitialConditions.Add(Equation3.Model.InitialConditions[0]);
			}
			(analyzers[2], solvers[2], nlAnalyzers[2]) = Equation3.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);


			for (int i = 0; i < analyzers.Length; i++)
            {
                analyzers[i].Initialize(true);
                if (analyzerStates[i] != null)
                {
                    analyzers[i].CurrentState = analyzerStates[i];
                }

				if (nlAnalyzerStates[i] != null)
				{
					nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
				}
			}
       }

       public void CreateModelFirstTime(IParentAnalyzer[] analyzers, ISolver[] solvers)
       {
            if (!(CurrentTimeStep == 0))
            {
				Equation1.UpdatePressureDivergenceDictionary(thetaTube, ParentSolvers[0], NLAnalyzers[0], model[0], Equation1.algebraicModel);
				Equation2.UpdatePressureDivergenceDictionary(thetaShell, ParentSolvers[1], NLAnalyzers[1], model[1], Equation2.algebraicModel);
				Equation3.UpdatePressureDivergenceDictionary(thetaWall, ParentSolvers[2], NLAnalyzers[2], model[2], Equation3.algebraicModel);
				
				Equation1.UpdateIndependedVariable(new Dictionary<int, double>[] { thetaWall });
				Equation2.UpdateIndependedVariable(new Dictionary<int, double>[] { thetaWall });
				Equation3.UpdateIndependedVariable(new Dictionary<int, double>[] { thetaTube, thetaShell });

				//Equation1.UpdateDependedVariable(thetaTube);
				//Equation2.UpdateDependedVariable(thetaShell);
				//Equation3.UpdateDependedVariable(thetaWall);
			}


            model = new Model[3];
            
           //Create Initial Model eq78 (fluid pressure)
            model[0] = Equation1.GetModel();
			Equation1.AddBoundaryConditions(model[0]);
			if (CurrentTimeStep == 0)
           {
				model[0].InitialConditions.Add(Equation1.Model.InitialConditions[0]);
           }
           (analyzers[0], solvers[0], nlAnalyzers[0]) = Equation1.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

           //Create model for eq9 (hyperelastic material)
           model[1] = Equation2.GetModel();
			Equation2.AddBoundaryConditions(model[1]);
			if (CurrentTimeStep == 0)
			{
				model[1].InitialConditions.Add(Equation2.Model.InitialConditions[0]);
			}
			(analyzers[1], solvers[1], nlAnalyzers[1]) = Equation2.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

			model[2] = Equation3.GetModel();
			Equation3.AddBoundaryConditions(model[2]);
			if (CurrentTimeStep == 0)
			{
				model[2].InitialConditions.Add(Equation3.Model.InitialConditions[0]);
			}
			(analyzers[2], solvers[2], nlAnalyzers[2]) = Equation3.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

			for (int i = 0; i < analyzers.Length; i++)
           {
               analyzers[i].Initialize(true);
               if (analyzerStates[i] != null)
               {
                   analyzers[i].CurrentState = analyzerStates[i];
               }

				if (nlAnalyzerStates[i] != null)
				{
					nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
				}
			}
       }

       public void SaveStateFromElements(List<IGlobalVector> resultsTube, List<IGlobalVector> resultsShell, List<IGlobalVector> resultsWall)
       {
			resultsTube.Add(this.parentAnalyzers[0].CurrentAnalysisResult);
			resultsShell.Add(this.parentAnalyzers[1].CurrentAnalysisResult);
			resultsWall.Add(this.parentAnalyzers[2].CurrentAnalysisResult);
			//Equation1.SaveStateFromElements(model[0]);
   //        Equation2.SaveStateFromElements(model[1]);
		 //  Equation3.SaveStateFromElements(model[2]);
		}
   }
}
