import java.io.*;
import java.util.Scanner;
import java.util.Arrays;

import ij.*;
import ij.io.*;
import ij.ImagePlus;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;

public class K_Nearest implements PlugInFilter {
    ImagePlus reference;        // Reference image
    int k;                      // Number of nearest neighbors
    int dLevel;                  // Wavelet decoposition level

	/*********************************************************/

    public int setup(String arg, ImagePlus imp) {
        reference = imp;
        ImageConverter ic = new ImageConverter(imp);
        ic.convertToGray8();
        return DOES_ALL;
    }
	
	/*********************************************************/

    public void run(ImageProcessor img) {

        GenericDialog gd = new GenericDialog("k-nearest neighbor search", IJ.getInstance());
        gd.addNumericField("Number of nearest neighbors (K):", 3, 0);
        gd.addNumericField("Wavelet decomposition level:", 1, 0);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        k = (int) gd.getNextNumber();
        dLevel = (int) gd.getNextNumber();

        SaveDialog sd = new SaveDialog("Open search folder...", "any file (required)", "");
        if (sd.getFileName() == null) return;
        String dir = sd.getDirectory();

        String[] fileNames = search(dir);

        ImageAccess image = new ImageAccess(img);
        double[][] dt = readDataSet(dLevel);
		

		double[][] featAux = new double[1][dLevel * 3 * 2 + 2];

		featAux[0] = Wavelet_Haar.waveletHaar(image, dLevel);
		
		double[] featureVector = Wavelet_Haar.normalize(featAux, false)[0];
				
		KNN(featureVector, dt, fileNames, k, dir);

    }
	
	/*********************************************************/
	
	
	public void KNN(double[] featureVector, double[][] dt, String[] fNames, int K, String dir) {
		
		double[] dist = getDistances(featureVector, dt);
		double[] distOrdered = dist;
		int n = dist.length;
		
		int[] indexOrdered = new int[n];
		double aux = 0;
		int iAux = 0;
		
		//initialize index array
		for(int i = 0; i < n; i++){
			indexOrdered[i] = i;			
		}
		
		//sort dist and indexes
		for(int i = 0; i < dist.length; i++){
			for(int j = i+1; j < dist.length; j++){
				if(distOrdered[j] < distOrdered[i]){
					aux = distOrdered[i];
					distOrdered[i] = distOrdered[j];
					distOrdered[j] = aux;
					
					iAux = indexOrdered[i];
					indexOrdered[i] = indexOrdered[j];
					indexOrdered[j] = iAux;		
				}
			}
		}
		
		int start = 0;
				
		double[] NearNeighDist  = Arrays.copyOfRange(distOrdered , start, K + start);
		int[] 	 NearNeighIndex = Arrays.copyOfRange(indexOrdered, start, K + start);
		String[] fNamesOrdered  = new String[K];
		
		
		for(int j = 0; j < NearNeighIndex.length; j++){
			int i = NearNeighIndex[j];
			fNamesOrdered[j] = fNames[i];
		}
				
		
		
		String fV = "";
		String fName;

		for(int j = 0; j < fNamesOrdered.length; j++) {			
			int i = NearNeighIndex[j];
			fName = fNamesOrdered[j];
			ImagePlus image = new Opener().openImage(dir, fName + ".gif");
			ImageAccess img = new ImageAccess(image.getProcessor());
			
			for(int k = 0; k < dt[0].length; k++){ fV += dt[i][k] + " "; }
			
			img.show(i + " - " + fName + ": " + fV + " | " + NearNeighDist[j]);	
			fV = "";
		}
		
	}
	
	/*********************************************************/
	
	public double[] getDistances(double[] featureVector, double[][] dt) {
		int ny = dt.length;
		int nx = dt[0].length;
		double[] distances = new double[ny];
		
		for(int i = 0; i < ny; i++){
			distances[i] = euclideanDistance(featureVector, dt[i]);
			
		}
		
		return distances;
	}
	
	/*********************************************************/
	
	public double euclideanDistance(double[] a1, double[] a2){
		int n = a1.length;
		
		double sum = 0;
		for (int i = 0; i < n; i++){
			sum += Math.pow(a1[i] - a2[i], 2);
		}
		return Math.sqrt(sum);		
	}
	
	/*********************************************************/

	public double[][] readDataSet(int dLevel) {
        try {
            
			//File dtFile = new File("//home//nfs//11369684//ImageJ//ImageJ//plugins//WaveletHaar//Dataset.txt");
			File dtFile = new File("C://Code//PI//ImageJ//plugins//WaveletHaar//Dataset.txt");
            
            
            Scanner sc = new Scanner(dtFile);
            int nx = dLevel * 3 * 2 + 2;
            int ny = 0;
			
            while(sc.hasNextLine()){
                ny++;
				sc.nextLine();
			}
			sc.close();
			
            // read in the data			
			double[][] dt = new double[ny][nx];
			String[] line;
			
            sc = new Scanner(dtFile);
            for(int u = 0; u < ny; u++) {
                if(sc.hasNextLine()){
					line = sc.nextLine().split(" ");
					for(int v = 0; v < nx; v++){
						dt[u][v] = Double.parseDouble(line[v]);
					}
                }
            }
			sc.close();
			
            return dt;

          } catch (Exception e) {
            e.printStackTrace();
          }
		  
		  double[][] error = {{0,0},
							  {0,0}};
		  
		  return error;
        
    }
	
	/*********************************************************/

    public String[] search(String dir) {     
        
        if (!dir.endsWith(File.separator))
            dir += File.separator;
        String[] list = new File(dir).list();  /* lista de arquivos */
        if (list == null) return list;

        double[][] featureMatrix = new double[list.length][dLevel * 3 * 2 + 2];
        ImageAccess input;
        

        for (int i=0; i<list.length; i++) {
            IJ.showStatus(i + "/" + list.length + ": " + list[i]);   /* mostra na interface */
            IJ.showProgress((double)i / list.length);  /* barra de progresso */
            File f = new File(dir + list[i]);

            if (!f.isDirectory()) {
                ImagePlus image = new Opener().openImage(dir, list[i]); /* abre imagem */
                if (image != null) {
                    input = new ImageAccess(image.getProcessor());
                    featureMatrix[i] = Wavelet_Haar.waveletHaar(input, dLevel);
					
					list[i] = list[i].split("\\.")[0];
                }
            }
        }

        Wavelet_Haar.createDataset(featureMatrix, list);

        IJ.showProgress(1.0);
        IJ.showStatus("");
				
		return list;

    }
}