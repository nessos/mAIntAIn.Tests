namespace MGroup.FEM.Thermal.Tests.ElectricalEquipment
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;

    using Xunit;

    public class ACGeneratorTest
	{
		[Fact]
		public void AlternatingCurrentGenerator()
		{
			double angleStep = 18;
			//var tuple = CalculateMechanicalInputValues(angleStep);
			//var rotation = tuple.Item1;
			//var angularVelocity = tuple.Item2;

			MechanicalInput mechanicalInput = new MechanicalInput(/*angleStep, rotation, angularVelocity*/);

			int numberOfTurns = 1;
			double magneticField = 2 * Math.Pow(10, -6);
			double area = 1.5 * 1.5;
			double resistance = 1;
			double angleDifference = Math.PI / 2;
			ACGenerator generator = new ACGenerator(numberOfTurns, magneticField, area, angleDifference, resistance, mechanicalInput);

			var dictionaryEMF = generator.CalculateEMFValuesByStep(angleStep/*, constant*/);
			var dictionaryFlux = generator.CalculateMagneticFluxValuesByStep(angleStep/*, constant*/);

		}

		public static Tuple<List<double>, List<double>> CalculateMechanicalInputValues(double angleStep)
		{
			double constant = 2; // custom angular velocity - ω = 2rad/s
			List<double> rotation = new List<double>();
			List<double> angularVelocity = new List<double>();
			double currentAngle = 0;
			while (currentAngle < 360)
			{
				rotation.Add(currentAngle);
				angularVelocity.Add(constant);
				currentAngle += angleStep;
			}
			return Tuple.Create(rotation, angularVelocity);
		}
	}
}
