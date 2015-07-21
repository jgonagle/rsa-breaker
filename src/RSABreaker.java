import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import Jama.Matrix;
import java.math.BigDecimal;

public class RSABreaker 
{
	//ordered collection of the first 50 primes
	private static List<Integer> primes = Arrays.asList(
		  	 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 
			 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
			 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,											 
			 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
			 179, 181, 191, 193, 197, 199, 211, 223, 227, 229);
	
	private int modulus;	
	private int smoothness;
	private int base;
	
	private int lastBasePower;
	
	private ArrayList<SmoothVector> smoothVectors;
	private SmoothVector targetSmoothVector;
	
	private double[][] linIndMatrix;
	private	double[][] targetVector;
	
	public RSABreaker(int modulus, int smoothness)
	{
		if ((modulus >= 1) && (smoothness >= 1))
		{
			System.out.println("Modulus: " + modulus);
			System.out.println("Smoothness: " + smoothness);
			
			this.modulus = modulus;		
			this.smoothness = smoothness;
			lastBasePower = 1;
			
			linIndMatrix = new double[smoothness][smoothness];
			targetVector = new double[smoothness][1];
			
			SmoothVector.smoothness = smoothness;
			
			smoothVectors = new ArrayList<SmoothVector>();
			
			Random rng = new Random();
			
			base = rng.nextInt(modulus - 1) + 1;
			
			while (!(coprimeToModulus(base) && isBRough(base)))
			{
				base = rng.nextInt(modulus - 1) + 1;
			}
			
			System.out.println("Base: " + base + "\n");
			
			findLinearlyIndependentVectors();
			
			Matrix A = new Matrix(linIndMatrix);
		    Matrix B = new Matrix(targetVector);
		    double[][] X = A.solve(B).getArray();
		    
		    double[] coefficientVector = new double[smoothness];
		    
		    for (int i = 0; i < smoothness; i++)
		    {
		    	coefficientVector[i] = X[i][0];
		    }
		    
		    System.out.println("\nCoefficient Matrix: ");
		    printVector(coefficientVector);
		    
		    
		    double numerator = 0 - targetSmoothVector.getBasePower();
		    double denominator = -1;
		    
		    for (int i = 0; i < smoothness; i++)
		    {
		    	numerator += X[i][0] * smoothVectors.get(i).getBasePower();
		    	denominator += X[i][0];
		    }
		    
		    double result = numerator / denominator;
		    
			System.out.println("\nOrder of " + base + " (mod " + modulus + ") = " + 
							   numerator + "\\" + denominator + " = " + result);
		}
	}

	private void findLinearlyIndependentVectors() 
	{
		int inverse;
		double[] newSmoothPowerVector;
		
		for (int i = 1; i < modulus; i++)
		{
			lastBasePower = multiply(lastBasePower, base, modulus);
			inverse = findInverse(lastBasePower, modulus);
			
			if ((newSmoothPowerVector = getSmoothVector(inverse)) != null)
			{
				if (areLinearIndependentVectors(smoothVectors, newSmoothPowerVector))
				{
					smoothVectors.add(new SmoothVector(i, inverse, newSmoothPowerVector));
					
					if (smoothVectors.size() == smoothness)
					{
						double[] targetSmoothPowerVector;
						
						for (int j = i + 1; j < modulus; j++)
						{
							lastBasePower = multiply(lastBasePower, base, modulus);
							inverse = findInverse(lastBasePower, modulus);
						
							if ((targetSmoothPowerVector = getSmoothVector(inverse)) != null)
							{
								targetSmoothVector = new SmoothVector(j, inverse, targetSmoothPowerVector);
								
								for (int k = 0; k < smoothness; k++)
								{
									targetVector[k][0] = targetSmoothPowerVector[k];
								}
								
								break;
							}
						}
						
						break;
					}
				}
			}
		}
		
		if (smoothVectors.size() != smoothness)
		{
			System.out.println("No linearly independent set of vectors found for " + smoothness + 
							   " smoothness on base " + base + ".  Choose another base");
			System.exit(0);
		}
		else
		{
			System.out.println("\n" + smoothness + " linearly independent smooth power vectors found!");
			
			for (int j = 0; j < smoothVectors.size(); j++)
			{
				double[] smoothPowerVector = smoothVectors.get(j).getSmoothPowerVector();
				
				System.out.print("vector_" + j + "(" + 
								 smoothVectors.get(j).getBasePower() + ") : ");
				printVector(smoothPowerVector);
			}
			
			System.out.println("\nTarget Vector (" +
							   targetSmoothVector.getBasePower() + ") : ");
			printVector(targetSmoothVector.getSmoothPowerVector());
		}
	}

	//assumption is that smoothVectors is linearly independent
	//then, show that newSmoothPowerVector is not in the span of smoothVectors' smooth power vectors
	private boolean areLinearIndependentVectors(ArrayList<SmoothVector> smoothVectors, 
												double[] newSmoothPowerVector) 
	{
		LIMatrix smoothVectorMatrix = new LIMatrix(smoothVectors, newSmoothPowerVector);
		linIndMatrix = smoothVectorMatrix.getMatrixArray();
		
		return smoothVectorMatrix.columnsLinearlyIndependent();
	}

	private boolean isBRough(int num) 
	{
		for (int i = 0; i < smoothness; i++)
		{
			if ((num % primes.get(i)) == 0)
			{
				return false;
			}
		}
		
		return true;
	}
	
	//returns vector of smooth prime exponents of num if num is b-smooth
	//else returns null (doubles as a smoothness test)
	private double[] getSmoothVector(int num)
	{
		int remainder = num;
		double[] smoothVector = new double[smoothness];
		
		for (int i = 0; i < smoothness; i++)
		{
			while ((remainder % primes.get(i)) == 0)
			{
				smoothVector[i]++;
				remainder /= primes.get(i);
			}
		}
		
		if (remainder != 1)
		{
			return null;
		}
		else
		{
			return smoothVector;
		}
	}

	private static int exponent(int num, int power, int someModulus) 
	{
		int answer = 1;
		
		while (power > 0)
		{
			if ((power & 1) == 1)
			{
				answer = multiply(answer, num, someModulus);
			}
			
			power = power >> 1;
			num = square(num, someModulus);
		}
		
		return answer;
	}

	private boolean coprimeToModulus(int num) 
	{
		return areCoprime(num, modulus);
	}
	
	private static int[] findBezoutPair(int numOne, int numTwo, int someModulus)
	{
		if (areCoprime(numOne, numTwo))
		{
			int prevX = 1;
			int prevY = 0;
			int prevRemainder = numOne;
			
			int curX = 0;
			int curY = 1;
			int curRemainder = numTwo;
			
			int quotient = 0; 
			
			while (curRemainder > 0)
			{
				int tempX = curX;
				int tempY = curY;
				int tempRemainder = curRemainder;
				
				quotient = prevRemainder / curRemainder;
				curRemainder = prevRemainder - (quotient * curRemainder);
				
				curX = prevX - (quotient * curX);
				curY = prevY - (quotient * curY);
				
				prevX = tempX;
				prevY = tempY;
				prevRemainder = tempRemainder;
			}
			
			return (new int[]{makePositive(prevX, someModulus), 
							  makePositive(prevY, someModulus)});
		}
		else
		{
			System.out.println("Cannot find Bezout pair. " + numOne + " and " +
							   numTwo + " are not coprime");
			
			return null;
		}
	}
	
	private static int findInverse(int num, int someModulus)
	{
		if (num == 1)
		{
			return 1;
		}
		if (areCoprime(num, someModulus))
		{
			return (((someModulus * findInverse((someModulus % num), num)) + 1) / num);
		}
		else
		{
			System.out.println("Cannot invert " + num + " (mod " + someModulus +  ")");
			System.exit(0);
			return 0;
		}
	}
	
	/*
	private static int findInverse(int num, int someModulus)
	{
		if (areCoprime(num, someModulus))
		{
			return findBezoutPair(num, someModulus, someModulus)[0];
		}
		else
		{
			System.out.println("Cannot invert " + num + " (mod " + someModulus +  ")");
			return 0;
		}
	}
	*/
	
	private static boolean areCoprime(int numOne, int numTwo) 
	{
		int remainderOne = Math.min(numOne, numTwo);
		int remainderTwo = Math.max(numOne, numTwo);
		
		while (remainderOne != 0)
		{
			int temp = remainderOne;
			
			remainderOne = remainderTwo % remainderOne;
			remainderTwo = temp;
		}
		
		return (remainderTwo == 1);
	}

	private static int makePositive(int num, int someModulus)
	{
		return (((num % someModulus) + someModulus) % someModulus);
	}
	
	private static int add(int numOne, int numTwo, int someModulus)
	{
		return makePositive((numOne + numTwo), someModulus);
	}
	
	private static int subtract(int numOne, int numTwo, int someModulus)
	{
		return makePositive((numOne - numTwo), someModulus);
	}
	
	private static int multiply(int numOne, int numTwo, int someModulus)
	{
		return makePositive((numOne * numTwo), someModulus);
	}
	
	private static int divide(int numOne, int numTwo, int someModulus)
	{
		return multiply(numOne, findInverse(makePositive(numTwo, someModulus), someModulus), someModulus);
	}
	
	private static int square(int num, int someModulus)
	{
		return multiply(num, num, someModulus);
	}
	
	private void printVector(double[] vector)
	{
		System.out.print("[");
		
		for (int h = 0; h < vector.length; h++)
		{
			System.out.print(vector[h]);
			
			if ((h + 1) == vector.length)
			{
				System.out.println("]");
			}
			else
			{
				System.out.print(",");
			}
		}
	}

	public static void main(String[] args)
	{
		RSABreaker semiPrime = new RSABreaker(47053, 3);
	}
}