namespace MGroup.FEM.Thermal.Tests
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	public class Alternator
	{
		private readonly double frequency;
		private readonly double inertiaConstant;
		private readonly double[] inputPower;
		private readonly double[] current = new double[2];
		private readonly double busBarVoltage;
		private readonly double reactancePreFault;
		private readonly double reactanceDuringFault;
		private readonly double reactanceAfterFault;
		private readonly double[] transientInernalVoltage = new double[2];
		private readonly double maxPowerPreFault;
		private readonly double maxPowerDuringFault;
		private readonly double maxPowerAfterFault;

		public double Frequency { get => frequency; }
		public double InertiaConsant { get => inertiaConstant; }
		public double BusBarVoltage { get => busBarVoltage; }
		public double ReactancePreFault { get => reactancePreFault; }
		public double ReactanceDuringFault { get => reactanceDuringFault; }
		public double RreactanceAfterFault { get => reactanceAfterFault; }
		public double MaxPowerPreFault { get => maxPowerPreFault; }
		public double MaxPowerDuringFault { get => maxPowerDuringFault; }
		public double MaxPowerAfterFault { get => maxPowerAfterFault; }

		public Alternator(double frequency, double inertiaConstant, double[] inputPower, double busBarVoltage,
			double reactancePreFault, double reactanceDuringFault, double reactanceAfterFault)
		{
			this.frequency = frequency;
			this.inertiaConstant = inertiaConstant;
			this.inputPower = inputPower;
			this.busBarVoltage = busBarVoltage;
			this.reactancePreFault = reactancePreFault;
			this.reactanceDuringFault = reactanceDuringFault;
			this.reactanceAfterFault = reactanceAfterFault;

			current[0] = inputPower[0] / this.busBarVoltage;
			current[1] = inputPower[1] / this.busBarVoltage;

			transientInernalVoltage[0] = busBarVoltage + current[1] * reactancePreFault;
			transientInernalVoltage[1] = current[0] * reactancePreFault;

			maxPowerPreFault = Math.Sqrt(Math.Pow(transientInernalVoltage[0], 2) + Math.Pow(transientInernalVoltage[1], 2)) * Math.Abs(busBarVoltage) / reactancePreFault;
			maxPowerDuringFault = Math.Sqrt(Math.Pow(transientInernalVoltage[0], 2) + Math.Pow(transientInernalVoltage[1], 2)) * Math.Abs(busBarVoltage) / reactanceDuringFault;
			maxPowerAfterFault = Math.Sqrt(Math.Pow(transientInernalVoltage[0], 2) + Math.Pow(transientInernalVoltage[1], 2)) * Math.Abs(busBarVoltage) / reactanceAfterFault;
		}

		public double[] GetCurrent() => current;

		public double[] GetTransientInternalVoltage() => transientInernalVoltage;

		public double[] GetInputPower() => inputPower;

		public List<(double, double)> SolveSwingEquation(double timeStart, double timeEnd, double timeStep, double faultTimeStart, double faultTimeEnd)
		{
			List<(double, double)> angleInTime = new List<(double, double)>();
			double time = timeStart;
			double angle1 = Math.Atan(transientInernalVoltage[1] / transientInernalVoltage[0]) * 180 / Math.PI;
			double dAngle = 0;
			double angularMomentum = inertiaConstant / (180 * frequency);
			double dError = timeStep * Math.Pow(10, -3);
			while (time < (timeEnd + dError))
			{
				angle1 = angle1 + dAngle;

				double maxPower;
				if (time < faultTimeStart - dError)
				{
					maxPower = maxPowerPreFault;
				}
				else if ((time > faultTimeStart - dError) && (time < faultTimeStart + dError))
				{
					double power1 = inputPower[0] - (maxPowerPreFault * Math.Sin(angle1 * Math.PI / 180));
					double power2 = inputPower[0] - (maxPowerDuringFault * Math.Sin(angle1 * Math.PI / 180));
					double averagePower = (power1 + power2) / 2;
					dAngle += Math.Pow(timeStep, 2) * averagePower / angularMomentum;
					angleInTime.Add((angle1, time));
					time += timeStep;
					continue;
				}
				else if (time < faultTimeEnd - dError)
				{
					maxPower = maxPowerDuringFault;
				}
				else if ((time > faultTimeEnd - dError) && (time < faultTimeEnd + dError))
				{
					double power1 = inputPower[0] - (maxPowerDuringFault * Math.Sin(angle1 * Math.PI / 180));
					double power2 = inputPower[0] - (maxPowerAfterFault * Math.Sin(angle1 * Math.PI / 180));
					double averagePower = (power1 + power2) / 2;
					dAngle += Math.Pow(timeStep, 2) * averagePower / angularMomentum;
					angleInTime.Add((angle1, time));
					time += timeStep;
					continue;
				}
				else
				{
					maxPower = maxPowerAfterFault;
				}

				double electricalPower = maxPower * Math.Sin(angle1 * Math.PI / 180);
				double powerGain = inputPower[0] - electricalPower;
				dAngle += Math.Pow(timeStep, 2) * powerGain / angularMomentum;

				angleInTime.Add((angle1, time));
				time += timeStep;
			}
			return angleInTime;
		}
	}
}
