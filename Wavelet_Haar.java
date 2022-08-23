import java.io.File;
import java.io.FileWriter;
import java.text.Normalizer;
import java.util.Arrays;
import java.util.Collections;
import java.util.Vector;
import java.util.Locale;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import ij.*;


public class Wavelet_Haar {
    
    
    public static final int DECOMPOSITION_LIMIT = 10;
    
    public static double maxEntropy;
    public static double minEntropy;
    public static double maxEnergy;
    public static double minEnergy;

    public static double[] waveletHaar(ImageAccess img, int dLevel) {
        
        double[][] LL = img.getArrayPixels();

        return waveletHaar(LL, dLevel);

    }
	
	/*********************************************************/

    public static double[] waveletHaar(double[][] LL, int dLevel) {    
    
        int i = 0;
        int nFeatures = dLevel * 3 * 2 + 2;

        double[] featureVector = new double[nFeatures];
        double[] energyVector = new double[nFeatures/2];
        double[] entropyVector = new double[nFeatures/2];

        int ny = LL.length;
        int nx = LL[0].length;
        
        double L[][];
        double H[][];
        double LH[][];
        double HL[][];
        double HH[][];

        
        for(int k = 0; k < dLevel; k++){
            
            L = new double[ny][nx/2];
            H = new double[ny][nx/2];

            L = lowPass(LL, 0);
            H = highPass(LL, 0);

            LL = new double[ny/2][nx/2];
            LH = new double[ny/2][nx/2];
            HL = new double[ny/2][nx/2];
            HH = new double[ny/2][nx/2];

            LL = lowPass(L, 1);
            LH = highPass(L, 1);
            HL = lowPass(H, 1);
            HH = highPass(H, 1);

			entropyVector[i] = entropy(LH);
            energyVector[i]  = energy(LH);
            i++;
            entropyVector[i] = entropy(HL);
            energyVector[i]  = energy(HL);
            i++;
            entropyVector[i] = entropy(HH);
            energyVector[i]  = energy(HH);
            i++;

            nx = nx/2;
            ny = ny/2;                
        }

        entropyVector[0] = entropy(LL);
        energyVector[0]  = energy(LL);


        for(i = 0; i < nFeatures/2; i++){
            featureVector[2*i] 	 = entropyVector[i];
            featureVector[2*i+1] = energyVector[i];
        }

        return featureVector;
    }
	
	/*********************************************************/

    public static void createDataset(double[][] featureMatrix, String[] nameList) {
        try{            
            DecimalFormatSymbols fS = new DecimalFormatSymbols(Locale.US);
            DecimalFormat df = new DecimalFormat("0.000", fS);			
			fS.setGroupingSeparator('.');
			
			File dt = new File("Dataset.txt");
            dt.createNewFile();                

            //FileWriter FW = new FileWriter("//home//nfs//11369684//ImageJ//ImageJ//plugins//WaveletHaar//Dataset.txt");
			FileWriter FW = new FileWriter("C://Code//PI//ImageJ//plugins//WaveletHaar//Dataset.txt");

            featureMatrix = normalize(featureMatrix, true);
			
            for (int i = 0; i < featureMatrix.length; i++){
                for (int j = 0; j < featureMatrix[i].length; j++) {
                    if(j == featureMatrix[i].length-1) {
                        FW.write(df.format(featureMatrix[i][j]) + "\n");
                    } else {
                        FW.write(df.format(featureMatrix[i][j]) + " ");
                    }
                }
            }

            FW.close();

        } catch (Exception e) {
            System.err.println(e);
        }
    }
	
	/*********************************************************/

    public static double[][] lowPass(double[][] matrix, int axis){
        int n = 0;
        int ny = matrix.length;
        int nx = matrix[0].length;
        double lowPass[][] = new double[1][1];

        switch(axis){
            case 0: //x axis
                n = nx/2;
                lowPass = new double[ny][n];

                for (int u = 0; u < ny; u++){
                    for(int v = 0; v < n; v++){
                        lowPass[u][v] = ((matrix[u][v*2] + matrix[u][v*2+1])/2) * Math.sqrt(2);
                    }
                }
            break;

            case 1: //y axis
                n = ny/2;
                lowPass = new double[n][nx];

                for (int u = 0; u < n; u++){
                    for(int v = 0; v < nx; v++){
                        lowPass[u][v] = ((matrix[u*2][v] + matrix[u*2+1][v])/2) * Math.sqrt(2);
                    }
                }
            break;            
        }

        return shift(lowPass);        
    }

	/*********************************************************/

    public static double[][] highPass(double[][] matrix, int axis){
        int n = 0;
        int ny = matrix.length;
        int nx = matrix[0].length;
        double highPass[][] = new double[1][1];

        switch(axis){
            case 0: //x axis
                n = nx/2;
                highPass = new double[ny][n];

                for (int u = 0; u < ny; u++){
                    for(int v = 0; v < n; v++){
                        highPass[u][v] = ((matrix[u][v*2] - matrix[u][v*2+1])/2) * Math.sqrt(2);
                    }
                }
            break;

            case 1: //y axis
                n = ny/2;
                highPass = new double[n][nx];

                for (int u = 0; u < n; u++){
                    for(int v = 0; v < nx; v++){
                        highPass[u][v] = ((matrix[u*2][v] - matrix[u*2+1][v])/2) * Math.sqrt(2);
                    }
                }
            break;
        }

        return shift(highPass);        
    }
	
	/*********************************************************/

    public static double[][] shift (double[][] matrix) {
        double min = findMin(matrix);

        if (min > 0) min = 0;

        int ny = matrix.length;
        int nx = matrix[0].length;
        double[][] shiftedMatrix = new double[ny][nx];

        for (int u = 0; u < ny; u++){
            for(int v = 0; v < nx; v++){
                shiftedMatrix[u][v] = matrix[u][v] + Math.abs(min);
                if (shiftedMatrix[u][v] > 255){
                    shiftedMatrix[u][v] = 255;
                } else if (shiftedMatrix[u][v] < 0){
					shiftedMatrix[u][v] = 0;
                }
            }
        }		

        return shiftedMatrix;
    }
	
	/*********************************************************/

    private static double findMin (double[][] matrix){
        double ny = matrix.length;   
        double nx = matrix[0].length;
        double min = matrix[0][0];

        for (int u = 0; u < ny; u++){
            for(int v = 0; v < nx; v++){
                if(min > matrix[u][v]){
                    min = matrix[u][v];
                }
            }
        }
        return min;
    }
		
	/*********************************************************/

    private static double findMinVec (double[] vector){
        double n = vector.length;
        double min = vector[0];

        for(int i = 1; i < n; i++){
            if(min > vector[i]){
                min = vector[i];
            }
        }

        return min;
    }
	
	/*********************************************************/

    private static double findMax (double[][] matrix){
        double ny = matrix.length;   
        double nx = matrix[0].length;
        double max = matrix[0][0];

        for (int u = 0; u < ny; u++){
            for(int v = 0; v < nx; v++){
                if(max < matrix[u][v]){
                    max = matrix[u][v];
                }
            }
        }
        return max;
    }
	
	/*********************************************************/

    private static double findMaxVec (double[] vector){
        double n = vector.length;
        double max = vector[0];

        for(int i = 1; i < n; i++){
            if(max < vector[i]){
                max = vector[i];
            }
        }

        return max;
    }
    
	/*********************************************************/

    public static double entropy(double[][] matrix){
        int ny = matrix.length;
        int nx = matrix[0].length;

        double sum = 0;
        double log = 0;

        for (int u = 0; u < ny; u++){
            for(int v = 0; v < nx; v++){
                log = Math.log(matrix[u][v]);
                if(Double.isNaN(log)){
                    log = 0;
                }
                sum += matrix[u][v] * log;
            }
        }

        return sum;
    }
	
	/*********************************************************/

    public static double energy(double[][] matrix){
        int ny = matrix.length;
        int nx = matrix[0].length;

        double sum = 0;

        for (int u = 0; u < ny; u++){
            for(int v = 0; v < nx; v++){
                sum -= Math.pow(matrix[u][v], 2);                
            }
        }

        return sum;
    }

    
	/*********************************************************/

    public static double[][] normalize(double[][] featMatrix, boolean isDataset){
        int ny = featMatrix.length;
        int nx = featMatrix[0].length;

        double[][] entropyMatrix = new double[ny][nx/2];
        double[][] energyMatrix = new double[ny][nx/2];

        for(int i = 0; i < ny; i++){
            for(int j = 0; j < nx/2; j++){
                entropyMatrix[i][j] = featMatrix[i][j*2];
                energyMatrix[i][j] = featMatrix[i][j*2+1];
            }
        }

        if(isDataset){

            minEntropy = findMin(entropyMatrix);
            maxEntropy = findMax(entropyMatrix);

            minEnergy = findMin(energyMatrix);
            maxEnergy = findMax(energyMatrix);
        }



        for(int i = 0; i < ny; i++){
            for(int j = 0; j < nx/2; j++){
                entropyMatrix[i][j] = ((entropyMatrix[i][j] - minEntropy) / (maxEntropy - minEntropy));
                entropyMatrix[i][j] = (double) Math.round(entropyMatrix[i][j] * 1000) / 1000;

                energyMatrix[i][j] = ((energyMatrix[i][j] - minEnergy) / (maxEnergy - minEnergy));
                energyMatrix[i][j] = (double) Math.round(energyMatrix[i][j] * 1000) / 1000;

            }
        }

        for(int i = 0; i < ny; i++){
            for(int j = 0; j < nx/2; j++){
                featMatrix[i][j*2] = entropyMatrix[i][j];
                featMatrix[i][j*2+1] = energyMatrix[i][j];
            }
        }

        return featMatrix;
    }
	
	/*********************************************************/

    public static double[] normalize(double[] vector){
        int n = vector.length;
        double aux;
        double min = findMinVec(vector);
        double max = findMaxVec(vector);
        
        for (int i = 0; i < n; i++) {
            aux = vector[i];
            vector[i] = ((vector[i] - min) / (max - min));
			vector[i] = (double) Math.round(vector[i] * 1000) / 1000;
        }

        return vector;
    }
	
	/*********************************************************/

    public static double[] getColumn(double[][] array, int col){
        double[] column = new double[array.length];

        for(int i = 0; i < column.length; i++){
           column[i] = array[i][col];
        }
        return column;
    }

    /*********************************************************/

    public static void printMatrix(double matrix[][]){
        int ny = matrix.length;
        int nx = matrix[0].length;

        System.out.println();

        for (int u = 0; u < ny; u++){
            for(int v = 0; v < nx; v++){
                IJ.log(Double.toString(matrix[u][v]) + " ");                
            }
            IJ.log("\n");
        }
    }
	
	/*********************************************************/
	
    public static void printVector(double vector[]){
        for(int i = 0; i < vector.length; i++){
            IJ.log(Double.toString(vector[i]) + " ");
        }
        IJ.log("\n");
    }
}