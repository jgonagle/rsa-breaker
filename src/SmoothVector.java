public class SmoothVector 
{
	protected static int smoothness;
	
	private int basePower;
	private int congruence;
	private double[] powerVector;
	
	public SmoothVector(int basePower, int congruence, double[] smoothPowerVector)
	{
		if (smoothPowerVector.length == smoothness)
		{
			this.basePower = basePower;
			this.congruence = congruence;
			this.powerVector = smoothPowerVector;
		}
		else
		{
			System.out.println("Smooth power vector must be of dimension " + smoothness + "!");
			System.exit(0);
		}
	}
		
	public int getBasePower()
	{
		return basePower;
	}
	
	public int getCongruence()
	{
		return congruence;
	}
	
	public double[] getSmoothPowerVector()
	{
		return powerVector;
	}
}
