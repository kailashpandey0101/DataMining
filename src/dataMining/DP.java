package dataMining;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.Scanner;

import weka.core.Instances;
import weka.core.converters.ConverterUtils.DataSource;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import static java.lang.Math.log;

//import weka.attributeSelection.AttributeSelection;
import weka.attributeSelection.OneRAttributeEval;
import weka.attributeSelection.Ranker;
import weka.filters.Filter; 
import weka.filters.supervised.attribute.AttributeSelection;


public class DP {
	static int k;
	static String[][]geneData2=null;  //TODO: changed
	static Instances ini_instances;
	 static Double[] classData1 ;
	 static String[][]geneData1=null; //changed
	 static String filename;
		//function to calculate the highest inforamtion gain by considering each split in a given arrray.
	    public static double[] splitArray(Double[] gData, Double[] cData, double eAll){
	            double maxGain= 0, maxSplit=0, maxIndx=0; int []s=new int[3];
	            for (int i=0 ; i< gData.length-1; i++){
	                int s1=0, s2=0;
	                double splitVal= (gData[i]+gData[i+1]) /2;
	                for (Double geneData1 : gData) {
	                    if (geneData1 < splitVal) {
	                        s1++;
	                    } else {
	                        s2++;
	                    }
	                }
	               
	                double p1= 0,n1=0,p2=0,n2=0;
	                for (int a=0; a<s1;a++){
	                   if (cData[a]==1){
	                       p1++;
	                   }
	                   else
	                       n1++;
	                }
	                for (int b=s1; b<s2+s1;b++){
	                   if (cData[b]==1){
	                       p2++;
	                   }
	                   else
	                       n2++;
	                }
	                float pS1 = (float) (log(p1/s1)/log(2));
	                if (Double.isInfinite(pS1)) {
	                    pS1 = 0;
	                }
	                float nS1 = (float) (log(n1/s1)/log(2));
	                if (Double.isInfinite(nS1)) {
	                    nS1 = 0;
	                }
	                float pS2 = (float) (log(p2/s2)/log(2));
	                if (Double.isInfinite(pS2)) {
	                    pS2 = 0;
	                }
	                float nS2 = (float) (log(n2/s2)/log(2));
	                if (Double.isInfinite(nS2)) {
	                    nS2 = 0;
	                }
	                double entropyS1 =-(p1/s1)*(pS1)-(n1/s1)*(nS1);
	                double entropyS2 =-(p2/s2)*(pS2)-(n2/s2)*(nS2);
	                double f1 = (double)s1/(s1+s2);
	                double f2 = (double)s2/(s1+s2);
	                double inSplit =entropyS1*(f1)+ entropyS2*(f2);
	                double inGain =eAll - inSplit;
	                if (maxGain<inGain){
	                    maxGain = inGain;
	                    maxIndx = i;
	                    maxSplit= splitVal;
	                    s[0]=s1;s[1]=s2;s[2]=i;
	                }
	              
	            }
	         
	            double[] maxInfo = new double[3];
	            maxInfo[0] = maxGain;
	            maxInfo[1]= maxIndx;
	            maxInfo[2]= maxSplit;
	           
	            return maxInfo;
	            
	    }
	    
		public static double calculateEntropy(Double[] classData1){
	        float cntZero = 0;
	        float cntOne = 0;
	        for(int i = 0; i<classData1.length; i++){
	            if (classData1[i]== 0){
	                cntZero+=1; 
	            }
	            else
	                cntOne+=1;
	        }
	        double p=cntOne/(cntZero+cntOne);
	        double n=cntZero/(cntZero+cntOne);
	        float pS = (float) (Math.log(p)/Math.log(2));
	                if (Double.isInfinite(pS)) {
	                    pS = 0;
	                }
	        float nS = (float) (Math.log(n)/Math.log(2));
	                if (Double.isInfinite(nS)) {
	                    nS = 0;
	                }
	        double entropyAll1 =(double)(-p*(pS)-n*(nS));
	        return entropyAll1;
	    }
	public static double getMean(double[] binmean)
	{
		double sum=0,average=0;
		for(int i=0;i<binmean.length;i++)
		{
			sum=sum+binmean[i];
		}
		average=sum/binmean.length;
		return average;
	}
	
	public static double getVariance(double[] binvar)
	{
		double mean=getMean(binvar);
		double temp=0;
		for (int i=0;i<binvar.length;i++)
		{
			temp+=(binvar[i]-mean)*(binvar[i]-mean);
			
		}
		return temp/binvar.length;
	}
	
	public static boolean hasData(double num, double[] bin)
	{
		boolean check=false;
		for(int i=0;i<bin.length;i++)
		{
			if (num==bin[i])
			{
				check=true;
			}
		}
		return check;
	}
	public static void eDensityBins (double[]oBin,double[] bin, String geneNo,int count,BufferedWriter bw) throws IOException{
		
		double[] firstBin,secondBin,thirdBin;
		firstBin=Arrays.copyOfRange(bin, 0, bin.length/3);
		secondBin=Arrays.copyOfRange(bin, firstBin.length, 2*(bin.length/3)+1);
		thirdBin=Arrays.copyOfRange(bin, 2*(bin.length/3)+1, bin.length);
		double a1=(bin[firstBin.length-1]+secondBin[0])/2;
		double a2=(bin[secondBin.length-1]+thirdBin[0])/2;
		for (int i=0;i<62;i++)
		{
			if(hasData(oBin[i],firstBin))
			{
				geneData2[i][count]="a";
			}
			else if(hasData(oBin[i],secondBin))
			{
				geneData2[i][count]="b";
			}
			else
			{
				geneData2[i][count]="c";
			}
		}
		
		bw.write(geneNo+" :");
		bw.write("variance :"+getVariance(firstBin)+", (- ,"+a1+"]"+" Count ="+firstBin.length+" ");
		bw.write("variance :"+getVariance(secondBin)+", ("+a1+","+a2+"]"+" Count ="+secondBin.length+" ");
		bw.write("variance :"+getVariance(thirdBin)+", ("+a2+",+]"+" Count ="+thirdBin.length);
		bw.newLine();
	}
	public static double Correlation(double[] xs, double[] ys) {
	    //TODO: check here that arrays are not null, of the same length etc

	    double sx = 0.0;
	    double sy = 0.0;
	    double sxx = 0.0;
	    double syy = 0.0;
	    double sxy = 0.0;

	    int n = (int)xs.length;

	    for(int i = 0; i < n; ++i) {
	      double x = xs[i];
	      double y = ys[i];

	      sx += x;
	      sy += y;
	      sxx += x * x;
	      syy += y * y;
	      sxy += x * y;
	    }

	    // covariation
	    double cov = sxy / n - sx * sy / n / n;
	    // standard error of x
	    double sigmax = Math.sqrt(sxx / n -  sx * sx / n / n);
	    // standard error of y
	    double sigmay = Math.sqrt(syy / n -  sy * sy / n / n);

	    // correlation is just a normalized covariation
	    return cov / sigmax / sigmay;
	  }
	private static void getBins(String geneName, int index,double topGain,int count,BufferedWriter bw) throws IOException {

		int gName = Integer.parseInt(stripNonDigits(geneName));
		double[] firstBin,secondBin;
		double a1;
		double geneData[]=new double[ini_instances.numInstances()];
		double unsortedGeneData[]=new double[ini_instances.numInstances()];

		for (int i=0; i<ini_instances.numInstances(); i++){
            geneData[i] = ini_instances.instance(i).value(gName);
            unsortedGeneData[i]= ini_instances.instance(i).value(gName);
        }
		Arrays.sort(geneData);
		firstBin=Arrays.copyOfRange(geneData, 0,index+1);
		secondBin=Arrays.copyOfRange(geneData, index+1,geneData.length);
		a1=round((geneData[firstBin.length-1]+secondBin[0])/2,4);
		bw.write(geneName+":"+"\tInfo Gain:  "+topGain+";\tBins: "+"(-, "+a1+" )\t"+(firstBin.length)+"\t\t("+a1 +",+)\t\t"+(secondBin.length));
		bw.newLine();// write to entropybins
		for(int i=0;i<62;i++)
		{
			if (hasData(unsortedGeneData[i],firstBin))
			{
				geneData1[i][count]="a";
			}
			else
			{
				geneData1[i][count]="b";
			}
		}
		for (int i=0;i<62;i++){
			if(classData1[i]==1)
			{
				geneData1[i][200]="positive";  //TODO: changed
			}
			else
			{
				geneData1[i][200]="negative";
			}
		}

		
	}
	public static String stripNonDigits(
            final CharSequence input ){
    final StringBuilder sb = new StringBuilder(
            input.length() );
    for(int i = 0; i < input.length(); i++){
        final char c = input.charAt(i);
        if(c > 47 && c < 58){
            sb.append(c);
        }
    }
    return sb.toString();
}

	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();

	    long factor = (long) Math.pow(10, places);
	    value = value * factor;
	    long tmp = Math.round(value);
	    return (double) tmp / factor;
	}
	
	
	public static void main(String args[]) throws Throwable{
		try
		{
			
			Scanner sc=new Scanner (System.in);
			System.out.println("Input the name of the file");
			filename=sc.next();
			System.out.println("Enter the value of k");
			k=Integer.parseInt(sc.next());
			geneData2=new String[62][k+1];
			geneData1=new String[62][k+1];
			
			File file4=new File(filename);
			Scanner inputFile =new Scanner(file4);
			File file2=new File(filename+".csv");
		
			file2.createNewFile();
			FileWriter writer=new FileWriter(file2);
			//CSVWriter writer=new CSVWriter(file2);
			while (inputFile.hasNext())
			{
				String csv=inputFile.nextLine();
				writer.append(csv);
				writer.append("\n");
				writer.flush();
			}
			DataSource fdata=new DataSource (filename+".csv");
			Instances _instances1=fdata.getDataSet();
			//Adding genes names to the given data
			String header ="G1";
			for(int i=2;i<_instances1.numAttributes();i++)
			{
				header=header+", G"+i;
			}
			header=header+", class";
			RandomAccessFile file1 = new RandomAccessFile(filename+".csv", "rws");
	        byte[] text = new byte[(int) file1.length()];
	        file1.readFully(text);
	        file1.seek(0);
	        file1.writeBytes(header);
	        file1.writeBytes("\n");
	        file1.write(text);
	        file1.close();
	        //inputFile.close();
	        //writer.close();
	        //sc.close();
	        
		}
		catch (Exception e){
			
		}
//		--------------Task 1 Starts--------------------
		System.out.println("-----------------------TASK 1 STARTS--------------------");
		//added from TASK 2
		Map<String,Double> indexSplit=new HashMap<String,Double>();
		Map<String,Double> geneCor= new HashMap<String,Double>();
		BufferedWriter bw_edensitybin = new BufferedWriter(new FileWriter("edensitybins.txt"));
		BufferedWriter bw_edensitydata = new BufferedWriter(new FileWriter("edensitydata.txt"));
		BufferedWriter bw_entropybin = new BufferedWriter(new FileWriter("entropybins.txt"));
		BufferedWriter bw_entropydata = new BufferedWriter(new FileWriter("entropydata.txt"));
		BufferedWriter bw_cor = new BufferedWriter(new FileWriter("correlatedgenes.txt"));
		
		DataSource frData;  
        frData = new DataSource( filename+".csv" );     //TODO: name to be changed
        ini_instances = frData.getDataSet();
        
        double[]iGain=new double[ini_instances.numAttributes()];
        
        ini_instances.setClassIndex( ini_instances.numAttributes() - 1 );
        Double[] geneData = new Double[ini_instances.numInstances()];
        Double[] classData = new Double[ini_instances.numInstances()];
        classData1 = new Double[ini_instances.numInstances()];
        // added from task 2
		File file=new File(filename+".csv");   //TODO: need to be changed
		double[][] myStringArray = new double [62][k];
		double[][] geneArray=new double[ini_instances.numInstances()][ini_instances.numAttributes()-1];
		Instances ini_instances;  
      
        ini_instances = frData.getDataSet();
		Scanner scan=new Scanner(file);
		int row=0;
		scan.nextLine();
		while (scan.hasNext()) {
            String csv = scan.nextLine();
            String[] spt=csv.split(",");
            for (int i=0;i<=k-1;i++)
            {
            	myStringArray[row][i]=Double.parseDouble(spt[i]);
            }
            row+=1;
        }
		for(int j = 0;j<ini_instances.numInstances(); j++ ){
            classData[j] = ini_instances.instance(j).value(ini_instances.numAttributes()-1);
            
        }
		double[] unSortGenes=new double[62];
		double[] sortGenes=new double[62];
		for(int i=0;i<k;i++)
		{
			for (int j=0;j<62;j++)
			{
				unSortGenes[j]=myStringArray[j][i];
				sortGenes[j]=myStringArray[j][i];
			}
			
			Arrays.sort(sortGenes);
			
			eDensityBins(unSortGenes,sortGenes,"G"+(i+1),i,bw_edensitybin);
		}
		bw_edensitybin.close();
		System.out.println("-----------------------edensitybin.txt file created--------------------");
		for (int i=0;i<62;i++){
			if (classData[i]==1){
				geneData2[i][k]="positive";
			}
			else
			{
				geneData2[i][k]="negative";
			}
		}

		for (int i=0;i<ini_instances.numInstances();i++)
		{
			for(int j=0;j<k+1;j++)
			{
				bw_edensitydata.write(geneData2[i][j]+" ");				// write to edensitydata.txt
			}
			bw_edensitydata.newLine();
		}
		bw_edensitydata.close();
		System.out.println("-----------------------edensitydata.txt file created--------------------");
//		---------------Task 1 Ends------------------
		System.out.println("--------------TASK 1 ENDS --------------------");
//		*********************Task 2 Starts*********************
		System.out.println("--------------TASK 2 STARTS --------------------");
		for(int j = 0;j<ini_instances.numInstances(); j++ ){
            classData[j] = ini_instances.instance(j).value(ini_instances.numAttributes()-1);
            classData1[j] = ini_instances.instance(j).value(ini_instances.numAttributes()-1);

        }
        double entropyAll = calculateEntropy(classData);//calling function to calcuate entropy of values for a single gene
        for (int j=0; j<ini_instances.numAttributes()-1; j++){
            for (int i=0; i<ini_instances.numInstances(); i++){
                geneData[i] = ini_instances.instance(i).value(j);
               
            }
            Map<Double, Double> myLocalMap = new HashMap<Double, Double>();
            for(int indx=0; indx < geneData.length; indx++){
                myLocalMap.put(Double.valueOf(geneData[indx]), Double.valueOf(classData[indx]));
             }
            Arrays.sort(geneData);

            for(int indx=0; indx < geneData.length; indx++){
                classData[indx] = myLocalMap.get(geneData[indx]).doubleValue();
            }
            double[] retSp = new double[3];
            retSp = splitArray(geneData, classData, entropyAll);//method called to find the highest information gain
            iGain[j]=round(retSp[0],6);									//Stores all the information gain for each gene
            
            indexSplit.put("G"+j+","+retSp[1], iGain[j]);					//stores info gain with respective split index
           
            
        }
		Map sortedMap = sortByValue(indexSplit);
		
		Set s = sortedMap.entrySet();
		
		String geneName="";
		int index;
		double topGain;
		int count=0;
		for(Object s1: s){
			String[] tokens = s1.toString().split("=");
			String [] spt=tokens[0].split(",");
			geneName=spt[0].trim();
			
			double ix=Double.parseDouble(spt[1].trim());
			index=(int)ix;
			topGain=Double.parseDouble(tokens[1].trim());
			if (count<k)
			{
				getBins(geneName,index,topGain,count,bw_entropybin);
				count++;
			}
			

		}
		bw_entropybin.close();
		System.out.println("-----------------------entropybin.txt file created--------------------");
		for (int i=0;i<ini_instances.numInstances();i++)
		{
			for(int j=0;j<201;j++)
			{
				bw_entropydata.write(geneData1[i][j]+" ");		//printing out the EntropyData.txt
			}
			bw_entropydata.newLine();
		}
		bw_entropydata.close();
		System.out.println("-----------------------Entropydata.txt file created--------------------");
		//**************Task 2 ends here******************
		
		System.out.println("--------------TASK 2 ENDS --------------------");
		//*************Task 3 starts**********************
		System.out.println("--------------TASK 3 STARTS --------------------");
		double[] arr1=new double[ini_instances.numInstances()];
		double[] arr2=new double[ini_instances.numInstances()];
		double cor;
		//System.out.println(ini_instances.numInstances());
		Scanner scan1=new Scanner(file);
		int row2=0;
		scan1.nextLine();
		while (scan1.hasNext()) {
            String csv = scan1.nextLine();
            String[] spt=csv.split(",");
            for (int i=0;i<=ini_instances.numAttributes()-2;i++)
            {
            	geneArray[row2][i]=Double.parseDouble(spt[i]);
            }
            row2+=1;
        }
//	
//		
		for(int i=0;i<ini_instances.numAttributes()-1;i++)
		{


			for(int j=i+1;j<ini_instances.numAttributes()-1;j++){
				for (int z=0;z<ini_instances.numInstances();z++){
					arr1[z]=geneArray[z][i];
					arr2[z]=geneArray[z][j];
				}
				cor=Correlation(arr1,arr2);
				
				geneCor.put("G"+i+":"+"G"+j, round(cor,5));
			}
		}
		Map sortCor=sortByValue(geneCor);
		Set cs=sortCor.entrySet();

		for(Object s1:cs)
		{
			String[] tokens = s1.toString().split("=");
			bw_cor.write(tokens[0]+"\t"+tokens[1]);				// write to correlatedgenes.txt
			bw_cor.newLine();
			
		}
		
		
		scan1.close();
		scan.close();
		bw_cor.close();
		System.out.println("-----------------------Correlatedgenes.txt File created--------------------");
		System.out.println("-----------------------TASK 3 ENDS--------------------");
	}

	public static Map sortByValue(Map unsortedMap) {
		Map sortedMap = new TreeMap(new ValueComparator(unsortedMap));
		sortedMap.putAll(unsortedMap);
		return sortedMap;
	}
	
}
class ValueComparator implements Comparator {
	Map map;
 
	public ValueComparator(Map map) {
		this.map = map;
	}
 
	public int compare(Object keyA, Object keyB) {
		Comparable valueA = (Comparable) map.get(keyA);
		Comparable valueB = (Comparable) map.get(keyB);
		return valueB.compareTo(valueA);
	}
}
