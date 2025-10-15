namespace MGroup.FEM.Thermal.Tests.ElectricalEquipment
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using Xunit;

	public class AlternatorTest
	{
		[Fact]
		public void Test()
		{
			double frequency = 50;
			double inertiaConstant = 5;
			double[] inputPower = new double[] { 0.8, 0.074 };
			double busBarVoltage = 1;
			double reactancePreFault = 0.65;
			double reactanceDuringFault = 1.8;
			double reactanceAfterFault = 0.8;
			var alternator = new Alternator(frequency, inertiaConstant, inputPower, busBarVoltage, reactancePreFault, reactanceDuringFault, reactanceAfterFault);

			double timeStart = 0;
			double timeEnd = 1;
			double timeStep = 0.05;

			double faultTimeStart = 0.0;
			double faultTimeEnd = 0.3;
			var results = alternator.SolveSwingEquation(timeStart, timeEnd, timeStep, faultTimeStart, faultTimeEnd);
			CheckResults(results);
		}

		public void CheckResults(List<(double, double)> results)
		{
			double[] calculatedResults = new double[results.Count];
			for (int i = 0; i < calculatedResults.Length; i++)
			{
				calculatedResults[i] = results[i].Item1;
			}
			double[] expectedResults = new double[]
				{26.3870, 27.5371, 30.9348, 36.4289, 43.7861, 52.7193, 62.9250, 72.4986, 79.3955, 83.4235,
					84.5135, 82.6524, 77.8641, 70.2416, 60.0252, 47.7078, 34.1222, 20.4446, 8.0683, -1.6318, -7.5445,};
			for (int i = 0; i < calculatedResults.Length; i++)
			{
				double error = Math.Abs((calculatedResults[i] - expectedResults[i]) / expectedResults[i]) * 100;
				if (error > Math.Pow( 10, -1))
				{
					Assert.True(false);
					break;
				}
			}
			Assert.True(true);
		}
	}
}
