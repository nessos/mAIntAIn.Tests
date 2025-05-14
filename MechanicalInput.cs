namespace MGroup.FEM.Thermal.Tests
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	public class MechanicalInput
	{
		//public readonly Dictionary<double, double> rotationDictionary = new Dictionary<double, double>();

		//public readonly Dictionary<double, double> angularVelocityDictionary = new Dictionary<double, double>();

		//public readonly List<double> time = new List<double>();

		public MechanicalInput(/*double angleStep, List<double> rotationPerStep, List<double> angularVelocityPerStep*/)
		{
			//double angle = 0;
			//int count = 0;
			//while (angle < 360)
			//{
			//	rotationDictionary.Add(angle, rotationPerStep[count]);
			//	angularVelocityDictionary.Add(angle, angularVelocityPerStep[count]);
			//	time.Add((rotationPerStep[count] * Math.PI / 180) / angularVelocityPerStep[count]);
			//	angle += angleStep;
			//	count++;
			//}
		}

		public double FindRotation(double time)
		{
			double rotation = FindVelocity(time) * time;
			return rotation;
			//double rotation1 = 0;
			//double deg1 = 0;
			//double rotation2 = 0;
			//double deg2 = 0;
			//foreach (double deg in rotationDictionary.Keys)
			//{
			//	if (deg == angle)
			//	{
			//		return rotationDictionary[deg] * Math.PI / 180;
			//	}
			//	else if (deg < angle)
			//	{
			//		rotation1 = rotationDictionary[deg] * Math.PI / 180;
			//		deg1 = deg;
			//	}
			//	else
			//	{
			//		rotation2 = rotationDictionary[deg] * Math.PI / 180;
			//		deg2 = deg;
			//		break;
			//	}
			//}
			////double rotation = (rotation1 + rotation2) / 2;
			//double rotation = ((rotation1 * (angle - deg2)) - (rotation2 * (angle - deg1))) / (deg1 - deg2);
			//return rotation;
		}

		public double FindVelocity(double time)
		{
			double velocity = 2;
			return velocity;
			//double velocity1 = 0;
			//double deg1 = 0;
			//double velocity2 = 0;
			//double deg2 = 0;
			//foreach (double deg in angularVelocityDictionary.Keys)
			//{
			//	if (deg == angle)
			//	{
			//		return angularVelocityDictionary[deg];
			//	}
			//	else if (deg < angle)
			//	{
			//		velocity1 = angularVelocityDictionary[deg];
			//		deg1 = deg;
			//	}
			//	else
			//	{
			//		velocity2 = angularVelocityDictionary[deg];
			//		deg2 = deg;
			//		break;
			//	}
			//}
			////double velocity = (velocity1 + velocity2) / 2;
			//double velocity = ((velocity1 * (angle - deg2)) - (velocity2 * (angle - deg1))) / (deg1 - deg2);
			//return velocity;
		}
	}
}
