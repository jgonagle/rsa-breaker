import java.util.ArrayList;


public class LIMatrix 
{
	private int rowCount;
	private int columnCount;	
	private double[][] matrix;
	private double[][] reducedMatrix;
	
	public LIMatrix(double[][] matrix)
	{
		rowCount = matrix.length;
		columnCount = matrix[0].length;
		this.matrix = matrix;
		
		reducedMatrix = new double[rowCount][columnCount];
		
		for (int row = 0; row < rowCount; row++)
		{
			for (int column = 0; column < columnCount; column++)
			{
				reducedMatrix[row][column] = matrix[row][column];
			}
		}
		
		reduceMatrix();
	}
	
	public LIMatrix(ArrayList<SmoothVector> smoothVectors, double[] newSmoothPowerVector)
	{
		rowCount = newSmoothPowerVector.length;
		columnCount = smoothVectors.size() + 1;
		
		if (rowCount >= columnCount)
		{
			matrix = new double[rowCount][columnCount];
			reducedMatrix = new double[rowCount][columnCount];
			
			double[] smoothPowerVector;
			
			for (int column = 0; column < columnCount; column++)
			{
				if (column < (columnCount - 1))
				{
					smoothPowerVector = smoothVectors.get(column).getSmoothPowerVector();
				}
				else
				{
					smoothPowerVector = newSmoothPowerVector;
				}
				
				for (int row = 0; row < rowCount; row++)
				{
					matrix[row][column] = smoothPowerVector[row];
					reducedMatrix[row][column] = smoothPowerVector[row];
				}
			}
			
			reduceMatrix();
		}
		else
		{
			System.out.println("Number of rows must be greater than number of columsn!");
			System.exit(0);
		}
	}

	private void reduceMatrix() 
	{
		double upperCoefficient, lowerCoefficient;
		
		for (int i = 0; i < columnCount - 1; i++)
		{
			upperCoefficient = reducedMatrix[i][i];
			
			for (int j = i + 1; j < rowCount; j++)
			{
				lowerCoefficient = reducedMatrix[j][i];
				
				for (int h = i; h < columnCount; h++)
				{
					reducedMatrix[j][h] = (upperCoefficient * reducedMatrix[j][h]) -
										  (lowerCoefficient * reducedMatrix[i][h]);
				}
			}
		}
	}
	
	public boolean columnsLinearlyIndependent()
	{
		boolean emptyRow;
		
		for (int i = 0; i < rowCount; i++)
		{
			emptyRow = true;
			
			for (int j = 0; j < columnCount; j++)
			{				
				if (reducedMatrix[i][j] != 0)
				{
					emptyRow = false;
					break;
				}
			}
			
			if (emptyRow)
			{
				return false;
			}
		}
		
		return true;
	}
	
	public void printMatrix()
	{
		print2DArray(matrix);
	}
	
	public void printReducedMatrix()
	{
		print2DArray(reducedMatrix);
	}
	
	private static void print2DArray(double[][] array)
	{
		for (int i = 0; i < array.length; i++)
		{
			for (int j = 0; j < array[0].length; j++)
			{
				System.out.print(array[i][j]);
				
				if (j == (array[0].length - 1))
				{
					System.out.println("");
				}
				else
				{
					System.out.print(",");
				}
			}
		}
		
		System.out.println("");
	}
	
	public double[][] getMatrixArray()
	{
		return matrix;
	}
}
