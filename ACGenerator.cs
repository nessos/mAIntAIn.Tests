namespace MGroup.FEM.Thermal.Tests
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	internal class ACGenerator
	{
		public double EMF { get => emf; }
		public double MagneticFlux { get => magnetixFlux; }
		public double Current { get => current; }

		// Number of turns of wire in the loop
		private readonly int numberOfTurns;

		// Magnetic field
		private readonly double magneticField;

		// Area of coil
		private readonly double area;

		// Electromagnetic force or emf
		private double emf;

		// Electric current
		private double current;

		// Electric resistance
		private double resistance;

		// Magnetic flux through circuit
		private double magnetixFlux;

		// Angle difference of the vector magnetic field and area normal vector with handle bar rotation
		private double angleDifference;

		// Mechanical input used to operate the generator
		private MechanicalInput mechanicalInput;
		public ACGenerator(int numberOfTurns, double magneticField, double area, double angleDifference, double resistance, MechanicalInput mechanicalInput)
		{
			this.numberOfTurns = numberOfTurns;
			this.magneticField = magneticField;
			this.area = area;
			this.angleDifference = angleDifference;
			this.resistance = resistance;
			this.mechanicalInput = mechanicalInput;
		}

		public void CalculateEMF(double angle)
		{
			emf = numberOfTurns * magneticField * area * mechanicalInput.FindVelocity(angle) * Math.Sin((angle * Math.PI / 180) + angleDifference);
		}

		public void CalculateMagneticFlux(double angle)
		{
			magnetixFlux = (numberOfTurns * magneticField * area) * Math.Cos((angle * Math.PI / 180) + angleDifference);
		}

		public void CalculateCurrent()
		{
			current = emf / resistance;
		}

		public Dictionary<double, double> CalculateEMFValuesByStep(double degStep)
		{
			Dictionary<double, double> EMFValuesPerAngle = new Dictionary<double, double>();
			double theta = 0;
			//foreach (var angle in mechanicalInput.angularVelocityDictionary.Keys)
			//{
			//	CalculateEMF(angle);
			//	EMFValuesPerAngle.Add(angle, EMF);
			//}
			//return EMFValuesPerAngle;

			while (theta < 360)
			{
				CalculateEMF(theta);
				EMFValuesPerAngle.Add(theta, EMF);
				theta += degStep;
			}
			return EMFValuesPerAngle;
		}

		public Dictionary<double, double> CalculateMagneticFluxValuesByStep(double degStep/*, double constant*/)
		{
			Dictionary<double, double> magneticFlxValuesPerAngle = new Dictionary<double, double>();
			//foreach (var angle in mechanicalInput.rotationDictionary.Keys)
			//{
			//	CalculateMagneticFlux(angle);
			//	magneticFlxValuesPerAngle.Add(angle, MagneticFlux);
			//}
			//return magneticFlxValuesPerAngle;


			double theta = 0;
			while (theta < 360)
			{
				CalculateMagneticFlux(theta);
				magneticFlxValuesPerAngle.Add(theta, MagneticFlux);
				theta += degStep;
			}
			return magneticFlxValuesPerAngle;
		}
	}
}
